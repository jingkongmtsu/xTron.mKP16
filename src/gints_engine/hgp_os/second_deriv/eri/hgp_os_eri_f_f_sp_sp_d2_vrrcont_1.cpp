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
#include "localmemscr.h"
using namespace localmemscr;

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
// BRA1 as redundant position, total RHS integrals evaluated as: 827952
// BRA2 as redundant position, total RHS integrals evaluated as: 827952
// KET1 as redundant position, total RHS integrals evaluated as: 704282
// KET2 as redundant position, total RHS integrals evaluated as: 704282
// the redundant position is: KET2
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
// BRA1  BRA2
// X  X
// X  Y
// X  Z
// Y  X
// Y  Y
// Y  Z
// Z  X
// Z  Y
// Z  Z
// ####
// BRA1  KET1
// X  X
// X  Y
// X  Z
// Y  X
// Y  Y
// Y  Z
// Z  X
// Z  Y
// Z  Z
// ####
// BRA2  BRA2
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####
// BRA2  KET1
// X  X
// X  Y
// X  Z
// Y  X
// Y  Y
// Y  Z
// Z  X
// Z  Y
// Z  Z
// ####
// KET1  KET1
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####

void hgp_os_eri_f_f_sp_sp_d2_vrrcont_1(const Double& ic2, const Double& jc2, const Double& jc2_1, const Double& jc2_2, const Double& jc2_3, const Double& alpha, const Double& beta, const Double& gamma, const Double& delta, const Double* SQ_ERI_G_S_S_S_vrr_array, const Double* SQ_ERI_F_S_S_S_vrr_array, const Double* SQ_ERI_G_S_P_S_vrr_array, const Double* SQ_ERI_F_S_P_S_vrr_array, const Double* SQ_ERI_G_S_S_P_vrr_array, Double* SQ_ERI_G_S_S_S_C3003, Double* SQ_ERI_F_S_S_S_C3003, Double* SQ_ERI_G_S_P_S_C1003003, Double* SQ_ERI_F_S_P_S_C1003003, Double* SQ_ERI_G_S_S_P_C1000003003, Double* SQ_ERI_G_S_S_S_C3003_a, Double* SQ_ERI_G_S_P_S_C1003003_a, Double* SQ_ERI_G_S_S_P_C1000003003_a, Double* SQ_ERI_G_S_P_S_C3003_c, Double* SQ_ERI_F_S_P_S_C3003_c, Double* SQ_ERI_G_S_S_S_C1003003, Double* SQ_ERI_F_S_S_S_C1003003, Double* SQ_ERI_G_S_S_P_C1001003003, Double* SQ_ERI_F_S_S_S_C3003_a, Double* SQ_ERI_F_S_P_S_C1003003_a, Double* SQ_ERI_G_S_P_S_C3003_ac, Double* SQ_ERI_G_S_S_S_C1003003_a, Double* SQ_ERI_G_S_S_P_C1001003003_a, Double* SQ_ERI_G_S_S_S_C3003_b, Double* SQ_ERI_F_S_S_S_C3003_b, Double* SQ_ERI_G_S_P_S_C1003003_b, Double* SQ_ERI_F_S_P_S_C1003003_b, Double* SQ_ERI_G_S_S_P_C1000003003_b, Double* SQ_ERI_G_S_S_S_C3003_c, Double* SQ_ERI_F_S_S_S_C3003_c, Double* SQ_ERI_G_S_P_S_C1003003_c, Double* SQ_ERI_F_S_P_S_C1003003_c, Double* SQ_ERI_G_S_S_P_C1000003003_c, Double* SQ_ERI_G_S_S_S_C3003_ab, Double* SQ_ERI_G_S_P_S_C1003003_ab, Double* SQ_ERI_G_S_S_P_C1000003003_ab, Double* SQ_ERI_G_S_P_S_C3003_bc, Double* SQ_ERI_F_S_P_S_C3003_bc, Double* SQ_ERI_G_S_S_S_C1003003_b, Double* SQ_ERI_F_S_S_S_C1003003_b, Double* SQ_ERI_G_S_S_P_C1001003003_b, Double* SQ_ERI_G_S_S_S_C3003_bb, Double* SQ_ERI_F_S_S_S_C3003_bb, Double* SQ_ERI_G_S_P_S_C1003003_bb, Double* SQ_ERI_F_S_P_S_C1003003_bb, Double* SQ_ERI_G_S_S_P_C1000003003_bb, Double* SQ_ERI_G_S_P_S_C1001003003, Double* SQ_ERI_F_S_P_S_C1001003003, Double* SQ_ERI_G_S_P_S_C1001003003_a, Double* SQ_ERI_G_S_P_S_C1000003003_c, Double* SQ_ERI_F_S_P_S_C1000003003_c, Double* SQ_ERI_F_S_P_S_C1001003003_a, Double* SQ_ERI_G_S_P_S_C1000003003_ac, Double* SQ_ERI_G_S_P_S_C1001003003_b, Double* SQ_ERI_F_S_P_S_C1001003003_b, Double* SQ_ERI_G_S_P_S_C1001003003_c, Double* SQ_ERI_F_S_P_S_C1001003003_c, Double* SQ_ERI_G_S_P_S_C1001003003_ab, Double* SQ_ERI_G_S_P_S_C1000003003_bc, Double* SQ_ERI_F_S_P_S_C1000003003_bc, Double* SQ_ERI_G_S_P_S_C1001003003_bb, Double* SQ_ERI_F_S_P_S_C1001003003_bb  ) 
{

  /************************************************************
   * transform the array form of integral into variables: SQ_ERI_G_S_S_S
   ************************************************************/
  Double I_ERI_G4x_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[0];
  Double I_ERI_G3xy_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[1];
  Double I_ERI_G3xz_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[2];
  Double I_ERI_G2x2y_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[3];
  Double I_ERI_G2xyz_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[4];
  Double I_ERI_G2x2z_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[5];
  Double I_ERI_Gx3y_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[6];
  Double I_ERI_Gx2yz_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[7];
  Double I_ERI_Gxy2z_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[8];
  Double I_ERI_Gx3z_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[9];
  Double I_ERI_G4y_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[10];
  Double I_ERI_G3yz_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[11];
  Double I_ERI_G2y2z_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[12];
  Double I_ERI_Gy3z_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[13];
  Double I_ERI_G4z_S_S_S_vrr = SQ_ERI_G_S_S_S_vrr_array[14];

  /************************************************************
   * transform the array form of integral into variables: SQ_ERI_F_S_S_S
   ************************************************************/
  Double I_ERI_F3x_S_S_S_vrr = SQ_ERI_F_S_S_S_vrr_array[0];
  Double I_ERI_F2xy_S_S_S_vrr = SQ_ERI_F_S_S_S_vrr_array[1];
  Double I_ERI_F2xz_S_S_S_vrr = SQ_ERI_F_S_S_S_vrr_array[2];
  Double I_ERI_Fx2y_S_S_S_vrr = SQ_ERI_F_S_S_S_vrr_array[3];
  Double I_ERI_Fxyz_S_S_S_vrr = SQ_ERI_F_S_S_S_vrr_array[4];
  Double I_ERI_Fx2z_S_S_S_vrr = SQ_ERI_F_S_S_S_vrr_array[5];
  Double I_ERI_F3y_S_S_S_vrr = SQ_ERI_F_S_S_S_vrr_array[6];
  Double I_ERI_F2yz_S_S_S_vrr = SQ_ERI_F_S_S_S_vrr_array[7];
  Double I_ERI_Fy2z_S_S_S_vrr = SQ_ERI_F_S_S_S_vrr_array[8];
  Double I_ERI_F3z_S_S_S_vrr = SQ_ERI_F_S_S_S_vrr_array[9];

  /************************************************************
   * transform the array form of integral into variables: SQ_ERI_G_S_P_S
   ************************************************************/
  Double I_ERI_G4x_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[0];
  Double I_ERI_G3xy_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[1];
  Double I_ERI_G3xz_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[2];
  Double I_ERI_G2x2y_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[3];
  Double I_ERI_G2xyz_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[4];
  Double I_ERI_G2x2z_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[5];
  Double I_ERI_Gx3y_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[6];
  Double I_ERI_Gx2yz_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[7];
  Double I_ERI_Gxy2z_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[8];
  Double I_ERI_Gx3z_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[9];
  Double I_ERI_G4y_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[10];
  Double I_ERI_G3yz_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[11];
  Double I_ERI_G2y2z_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[12];
  Double I_ERI_Gy3z_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[13];
  Double I_ERI_G4z_S_Px_S_vrr = SQ_ERI_G_S_P_S_vrr_array[14];
  Double I_ERI_G4x_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[15];
  Double I_ERI_G3xy_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[16];
  Double I_ERI_G3xz_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[17];
  Double I_ERI_G2x2y_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[18];
  Double I_ERI_G2xyz_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[19];
  Double I_ERI_G2x2z_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[20];
  Double I_ERI_Gx3y_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[21];
  Double I_ERI_Gx2yz_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[22];
  Double I_ERI_Gxy2z_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[23];
  Double I_ERI_Gx3z_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[24];
  Double I_ERI_G4y_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[25];
  Double I_ERI_G3yz_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[26];
  Double I_ERI_G2y2z_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[27];
  Double I_ERI_Gy3z_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[28];
  Double I_ERI_G4z_S_Py_S_vrr = SQ_ERI_G_S_P_S_vrr_array[29];
  Double I_ERI_G4x_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[30];
  Double I_ERI_G3xy_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[31];
  Double I_ERI_G3xz_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[32];
  Double I_ERI_G2x2y_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[33];
  Double I_ERI_G2xyz_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[34];
  Double I_ERI_G2x2z_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[35];
  Double I_ERI_Gx3y_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[36];
  Double I_ERI_Gx2yz_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[37];
  Double I_ERI_Gxy2z_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[38];
  Double I_ERI_Gx3z_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[39];
  Double I_ERI_G4y_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[40];
  Double I_ERI_G3yz_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[41];
  Double I_ERI_G2y2z_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[42];
  Double I_ERI_Gy3z_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[43];
  Double I_ERI_G4z_S_Pz_S_vrr = SQ_ERI_G_S_P_S_vrr_array[44];

  /************************************************************
   * transform the array form of integral into variables: SQ_ERI_F_S_P_S
   ************************************************************/
  Double I_ERI_F3x_S_Px_S_vrr = SQ_ERI_F_S_P_S_vrr_array[0];
  Double I_ERI_F2xy_S_Px_S_vrr = SQ_ERI_F_S_P_S_vrr_array[1];
  Double I_ERI_F2xz_S_Px_S_vrr = SQ_ERI_F_S_P_S_vrr_array[2];
  Double I_ERI_Fx2y_S_Px_S_vrr = SQ_ERI_F_S_P_S_vrr_array[3];
  Double I_ERI_Fxyz_S_Px_S_vrr = SQ_ERI_F_S_P_S_vrr_array[4];
  Double I_ERI_Fx2z_S_Px_S_vrr = SQ_ERI_F_S_P_S_vrr_array[5];
  Double I_ERI_F3y_S_Px_S_vrr = SQ_ERI_F_S_P_S_vrr_array[6];
  Double I_ERI_F2yz_S_Px_S_vrr = SQ_ERI_F_S_P_S_vrr_array[7];
  Double I_ERI_Fy2z_S_Px_S_vrr = SQ_ERI_F_S_P_S_vrr_array[8];
  Double I_ERI_F3z_S_Px_S_vrr = SQ_ERI_F_S_P_S_vrr_array[9];
  Double I_ERI_F3x_S_Py_S_vrr = SQ_ERI_F_S_P_S_vrr_array[10];
  Double I_ERI_F2xy_S_Py_S_vrr = SQ_ERI_F_S_P_S_vrr_array[11];
  Double I_ERI_F2xz_S_Py_S_vrr = SQ_ERI_F_S_P_S_vrr_array[12];
  Double I_ERI_Fx2y_S_Py_S_vrr = SQ_ERI_F_S_P_S_vrr_array[13];
  Double I_ERI_Fxyz_S_Py_S_vrr = SQ_ERI_F_S_P_S_vrr_array[14];
  Double I_ERI_Fx2z_S_Py_S_vrr = SQ_ERI_F_S_P_S_vrr_array[15];
  Double I_ERI_F3y_S_Py_S_vrr = SQ_ERI_F_S_P_S_vrr_array[16];
  Double I_ERI_F2yz_S_Py_S_vrr = SQ_ERI_F_S_P_S_vrr_array[17];
  Double I_ERI_Fy2z_S_Py_S_vrr = SQ_ERI_F_S_P_S_vrr_array[18];
  Double I_ERI_F3z_S_Py_S_vrr = SQ_ERI_F_S_P_S_vrr_array[19];
  Double I_ERI_F3x_S_Pz_S_vrr = SQ_ERI_F_S_P_S_vrr_array[20];
  Double I_ERI_F2xy_S_Pz_S_vrr = SQ_ERI_F_S_P_S_vrr_array[21];
  Double I_ERI_F2xz_S_Pz_S_vrr = SQ_ERI_F_S_P_S_vrr_array[22];
  Double I_ERI_Fx2y_S_Pz_S_vrr = SQ_ERI_F_S_P_S_vrr_array[23];
  Double I_ERI_Fxyz_S_Pz_S_vrr = SQ_ERI_F_S_P_S_vrr_array[24];
  Double I_ERI_Fx2z_S_Pz_S_vrr = SQ_ERI_F_S_P_S_vrr_array[25];
  Double I_ERI_F3y_S_Pz_S_vrr = SQ_ERI_F_S_P_S_vrr_array[26];
  Double I_ERI_F2yz_S_Pz_S_vrr = SQ_ERI_F_S_P_S_vrr_array[27];
  Double I_ERI_Fy2z_S_Pz_S_vrr = SQ_ERI_F_S_P_S_vrr_array[28];
  Double I_ERI_F3z_S_Pz_S_vrr = SQ_ERI_F_S_P_S_vrr_array[29];

  /************************************************************
   * transform the array form of integral into variables: SQ_ERI_G_S_S_P
   ************************************************************/
  Double I_ERI_G4x_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[0];
  Double I_ERI_G3xy_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[1];
  Double I_ERI_G3xz_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[2];
  Double I_ERI_G2x2y_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[3];
  Double I_ERI_G2xyz_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[4];
  Double I_ERI_G2x2z_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[5];
  Double I_ERI_Gx3y_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[6];
  Double I_ERI_Gx2yz_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[7];
  Double I_ERI_Gxy2z_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[8];
  Double I_ERI_Gx3z_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[9];
  Double I_ERI_G4y_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[10];
  Double I_ERI_G3yz_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[11];
  Double I_ERI_G2y2z_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[12];
  Double I_ERI_Gy3z_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[13];
  Double I_ERI_G4z_S_S_Px_vrr = SQ_ERI_G_S_S_P_vrr_array[14];
  Double I_ERI_G4x_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[15];
  Double I_ERI_G3xy_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[16];
  Double I_ERI_G3xz_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[17];
  Double I_ERI_G2x2y_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[18];
  Double I_ERI_G2xyz_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[19];
  Double I_ERI_G2x2z_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[20];
  Double I_ERI_Gx3y_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[21];
  Double I_ERI_Gx2yz_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[22];
  Double I_ERI_Gxy2z_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[23];
  Double I_ERI_Gx3z_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[24];
  Double I_ERI_G4y_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[25];
  Double I_ERI_G3yz_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[26];
  Double I_ERI_G2y2z_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[27];
  Double I_ERI_Gy3z_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[28];
  Double I_ERI_G4z_S_S_Py_vrr = SQ_ERI_G_S_S_P_vrr_array[29];
  Double I_ERI_G4x_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[30];
  Double I_ERI_G3xy_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[31];
  Double I_ERI_G3xz_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[32];
  Double I_ERI_G2x2y_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[33];
  Double I_ERI_G2xyz_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[34];
  Double I_ERI_G2x2z_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[35];
  Double I_ERI_Gx3y_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[36];
  Double I_ERI_Gx2yz_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[37];
  Double I_ERI_Gxy2z_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[38];
  Double I_ERI_Gx3z_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[39];
  Double I_ERI_G4y_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[40];
  Double I_ERI_G3yz_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[41];
  Double I_ERI_G2y2z_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[42];
  Double I_ERI_Gy3z_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[43];
  Double I_ERI_G4z_S_S_Pz_vrr = SQ_ERI_G_S_S_P_vrr_array[44];

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_S_C3003
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_S_C3003_coefs = ic2*jc2;
  SQ_ERI_G_S_S_S_C3003[0] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_G4x_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[1] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_G3xy_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[2] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_G3xz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[3] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_G2x2y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[4] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_G2xyz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[5] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_G2x2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[6] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_Gx3y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[7] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_Gx2yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[8] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_Gxy2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[9] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_Gx3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[10] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_G4y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[11] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_G3yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[12] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_G2y2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[13] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_Gy3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003[14] += SQ_ERI_G_S_S_S_C3003_coefs*I_ERI_G4z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3003
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_S_S_C3003_coefs = ic2*jc2;
  SQ_ERI_F_S_S_S_C3003[0] += SQ_ERI_F_S_S_S_C3003_coefs*I_ERI_F3x_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003[1] += SQ_ERI_F_S_S_S_C3003_coefs*I_ERI_F2xy_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003[2] += SQ_ERI_F_S_S_S_C3003_coefs*I_ERI_F2xz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003[3] += SQ_ERI_F_S_S_S_C3003_coefs*I_ERI_Fx2y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003[4] += SQ_ERI_F_S_S_S_C3003_coefs*I_ERI_Fxyz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003[5] += SQ_ERI_F_S_S_S_C3003_coefs*I_ERI_Fx2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003[6] += SQ_ERI_F_S_S_S_C3003_coefs*I_ERI_F3y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003[7] += SQ_ERI_F_S_S_S_C3003_coefs*I_ERI_F2yz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003[8] += SQ_ERI_F_S_S_S_C3003_coefs*I_ERI_Fy2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003[9] += SQ_ERI_F_S_S_S_C3003_coefs*I_ERI_F3z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1003003
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1003003_coefs = ic2*jc2_1;
  SQ_ERI_G_S_P_S_C1003003[0] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[1] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[2] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[3] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[4] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[5] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[6] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[7] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[8] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[9] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[10] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[11] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[12] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[13] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[14] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[15] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[16] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[17] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[18] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[19] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[20] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[21] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[22] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[23] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[24] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[25] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[26] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[27] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[28] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[29] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[30] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[31] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[32] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[33] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[34] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[35] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[36] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[37] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[38] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[39] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[40] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[41] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[42] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[43] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003[44] += SQ_ERI_G_S_P_S_C1003003_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1003003
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1003003_coefs = ic2*jc2_1;
  SQ_ERI_F_S_P_S_C1003003[0] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[1] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[2] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[3] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[4] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[5] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[6] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[7] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[8] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[9] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[10] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[11] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[12] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[13] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[14] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[15] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[16] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[17] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[18] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[19] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[20] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[21] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[22] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[23] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[24] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[25] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[26] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[27] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[28] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003[29] += SQ_ERI_F_S_P_S_C1003003_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_P_C1000003003
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_P_C1000003003_coefs = ic2*jc2_2;
  SQ_ERI_G_S_S_P_C1000003003[0] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G4x_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[1] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G3xy_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[2] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G3xz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[3] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2x2y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[4] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2xyz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[5] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2x2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[6] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gx3y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[7] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[8] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[9] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gx3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[10] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G4y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[11] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G3yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[12] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2y2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[13] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gy3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[14] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G4z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003[15] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G4x_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[16] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G3xy_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[17] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G3xz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[18] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2x2y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[19] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2xyz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[20] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2x2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[21] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gx3y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[22] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[23] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[24] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gx3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[25] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G4y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[26] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G3yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[27] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2y2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[28] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gy3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[29] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G4z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003[30] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G4x_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[31] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G3xy_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[32] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G3xz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[33] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[34] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[35] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[36] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[37] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[38] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[39] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[40] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G4y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[41] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G3yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[42] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[43] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003[44] += SQ_ERI_G_S_S_P_C1000003003_coefs*I_ERI_G4z_S_S_Pz_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_S_C3003_a
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_S_C3003_a_coefs = ic2*jc2*alpha;
  SQ_ERI_G_S_S_S_C3003_a[0] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_G4x_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[1] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_G3xy_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[2] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_G3xz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[3] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_G2x2y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[4] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_G2xyz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[5] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_G2x2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[6] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_Gx3y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[7] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_Gx2yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[8] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_Gxy2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[9] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_Gx3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[10] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_G4y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[11] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_G3yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[12] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_G2y2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[13] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_Gy3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_a[14] += SQ_ERI_G_S_S_S_C3003_a_coefs*I_ERI_G4z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1003003_a
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1003003_a_coefs = ic2*jc2_1*alpha;
  SQ_ERI_G_S_P_S_C1003003_a[0] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[1] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[2] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[3] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[4] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[5] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[6] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[7] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[8] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[9] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[10] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[11] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[12] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[13] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[14] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[15] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[16] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[17] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[18] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[19] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[20] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[21] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[22] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[23] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[24] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[25] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[26] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[27] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[28] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[29] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[30] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[31] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[32] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[33] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[34] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[35] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[36] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[37] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[38] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[39] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[40] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[41] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[42] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[43] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_a[44] += SQ_ERI_G_S_P_S_C1003003_a_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_P_C1000003003_a
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_P_C1000003003_a_coefs = ic2*jc2_2*alpha;
  SQ_ERI_G_S_S_P_C1000003003_a[0] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G4x_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[1] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G3xy_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[2] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G3xz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[3] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2x2y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[4] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2xyz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[5] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2x2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[6] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gx3y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[7] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[8] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[9] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gx3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[10] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G4y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[11] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G3yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[12] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2y2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[13] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gy3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[14] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G4z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[15] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G4x_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[16] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G3xy_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[17] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G3xz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[18] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2x2y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[19] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2xyz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[20] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2x2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[21] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gx3y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[22] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[23] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[24] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gx3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[25] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G4y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[26] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G3yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[27] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2y2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[28] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gy3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[29] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G4z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[30] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G4x_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[31] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G3xy_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[32] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G3xz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[33] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[34] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[35] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[36] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[37] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[38] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[39] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[40] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G4y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[41] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G3yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[42] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[43] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_a[44] += SQ_ERI_G_S_S_P_C1000003003_a_coefs*I_ERI_G4z_S_S_Pz_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C3003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C3003_c_coefs = ic2*jc2*gamma;
  SQ_ERI_G_S_P_S_C3003_c[0] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[1] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[2] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[3] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[4] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[5] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[6] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[7] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[8] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[9] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[10] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[11] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[12] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[13] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[14] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[15] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[16] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[17] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[18] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[19] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[20] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[21] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[22] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[23] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[24] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[25] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[26] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[27] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[28] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[29] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[30] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[31] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[32] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[33] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[34] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[35] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[36] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[37] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[38] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[39] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[40] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[41] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[42] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[43] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_c[44] += SQ_ERI_G_S_P_S_C3003_c_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C3003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C3003_c_coefs = ic2*jc2*gamma;
  SQ_ERI_F_S_P_S_C3003_c[0] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[1] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[2] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[3] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[4] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[5] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[6] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[7] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[8] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[9] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[10] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[11] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[12] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[13] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[14] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[15] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[16] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[17] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[18] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[19] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[20] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[21] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[22] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[23] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[24] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[25] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[26] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[27] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[28] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_c[29] += SQ_ERI_F_S_P_S_C3003_c_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_S_C1003003
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_S_C1003003_coefs = ic2*jc2_1;
  SQ_ERI_G_S_S_S_C1003003[0] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_G4x_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[1] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_G3xy_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[2] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_G3xz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[3] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_G2x2y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[4] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_G2xyz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[5] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_G2x2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[6] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_Gx3y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[7] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_Gx2yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[8] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_Gxy2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[9] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_Gx3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[10] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_G4y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[11] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_G3yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[12] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_G2y2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[13] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_Gy3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003[14] += SQ_ERI_G_S_S_S_C1003003_coefs*I_ERI_G4z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C1003003
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_S_S_C1003003_coefs = ic2*jc2_1;
  SQ_ERI_F_S_S_S_C1003003[0] += SQ_ERI_F_S_S_S_C1003003_coefs*I_ERI_F3x_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003[1] += SQ_ERI_F_S_S_S_C1003003_coefs*I_ERI_F2xy_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003[2] += SQ_ERI_F_S_S_S_C1003003_coefs*I_ERI_F2xz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003[3] += SQ_ERI_F_S_S_S_C1003003_coefs*I_ERI_Fx2y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003[4] += SQ_ERI_F_S_S_S_C1003003_coefs*I_ERI_Fxyz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003[5] += SQ_ERI_F_S_S_S_C1003003_coefs*I_ERI_Fx2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003[6] += SQ_ERI_F_S_S_S_C1003003_coefs*I_ERI_F3y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003[7] += SQ_ERI_F_S_S_S_C1003003_coefs*I_ERI_F2yz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003[8] += SQ_ERI_F_S_S_S_C1003003_coefs*I_ERI_Fy2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003[9] += SQ_ERI_F_S_S_S_C1003003_coefs*I_ERI_F3z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_P_C1001003003
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_P_C1001003003_coefs = ic2*jc2_3;
  SQ_ERI_G_S_S_P_C1001003003[0] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G4x_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[1] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G3xy_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[2] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G3xz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[3] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2x2y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[4] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2xyz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[5] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2x2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[6] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gx3y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[7] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[8] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[9] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gx3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[10] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G4y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[11] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G3yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[12] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2y2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[13] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gy3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[14] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G4z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003[15] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G4x_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[16] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G3xy_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[17] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G3xz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[18] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2x2y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[19] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2xyz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[20] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2x2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[21] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gx3y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[22] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[23] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[24] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gx3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[25] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G4y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[26] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G3yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[27] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2y2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[28] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gy3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[29] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G4z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003[30] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G4x_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[31] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G3xy_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[32] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G3xz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[33] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[34] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[35] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[36] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[37] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[38] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[39] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[40] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G4y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[41] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G3yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[42] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[43] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003[44] += SQ_ERI_G_S_S_P_C1001003003_coefs*I_ERI_G4z_S_S_Pz_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3003_a
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_S_S_C3003_a_coefs = ic2*jc2*alpha;
  SQ_ERI_F_S_S_S_C3003_a[0] += SQ_ERI_F_S_S_S_C3003_a_coefs*I_ERI_F3x_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_a[1] += SQ_ERI_F_S_S_S_C3003_a_coefs*I_ERI_F2xy_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_a[2] += SQ_ERI_F_S_S_S_C3003_a_coefs*I_ERI_F2xz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_a[3] += SQ_ERI_F_S_S_S_C3003_a_coefs*I_ERI_Fx2y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_a[4] += SQ_ERI_F_S_S_S_C3003_a_coefs*I_ERI_Fxyz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_a[5] += SQ_ERI_F_S_S_S_C3003_a_coefs*I_ERI_Fx2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_a[6] += SQ_ERI_F_S_S_S_C3003_a_coefs*I_ERI_F3y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_a[7] += SQ_ERI_F_S_S_S_C3003_a_coefs*I_ERI_F2yz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_a[8] += SQ_ERI_F_S_S_S_C3003_a_coefs*I_ERI_Fy2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_a[9] += SQ_ERI_F_S_S_S_C3003_a_coefs*I_ERI_F3z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1003003_a
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1003003_a_coefs = ic2*jc2_1*alpha;
  SQ_ERI_F_S_P_S_C1003003_a[0] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[1] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[2] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[3] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[4] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[5] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[6] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[7] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[8] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[9] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[10] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[11] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[12] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[13] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[14] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[15] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[16] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[17] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[18] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[19] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[20] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[21] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[22] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[23] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[24] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[25] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[26] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[27] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[28] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_a[29] += SQ_ERI_F_S_P_S_C1003003_a_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C3003_ac
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C3003_ac_coefs = ic2*jc2*alpha*gamma;
  SQ_ERI_G_S_P_S_C3003_ac[0] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[1] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[2] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[3] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[4] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[5] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[6] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[7] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[8] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[9] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[10] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[11] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[12] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[13] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[14] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[15] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[16] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[17] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[18] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[19] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[20] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[21] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[22] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[23] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[24] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[25] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[26] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[27] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[28] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[29] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[30] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[31] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[32] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[33] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[34] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[35] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[36] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[37] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[38] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[39] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[40] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[41] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[42] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[43] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_ac[44] += SQ_ERI_G_S_P_S_C3003_ac_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_S_C1003003_a
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_S_C1003003_a_coefs = ic2*jc2_1*alpha;
  SQ_ERI_G_S_S_S_C1003003_a[0] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_G4x_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[1] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_G3xy_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[2] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_G3xz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[3] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_G2x2y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[4] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_G2xyz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[5] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_G2x2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[6] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_Gx3y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[7] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_Gx2yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[8] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_Gxy2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[9] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_Gx3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[10] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_G4y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[11] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_G3yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[12] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_G2y2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[13] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_Gy3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_a[14] += SQ_ERI_G_S_S_S_C1003003_a_coefs*I_ERI_G4z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_P_C1001003003_a
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_P_C1001003003_a_coefs = ic2*jc2_3*alpha;
  SQ_ERI_G_S_S_P_C1001003003_a[0] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G4x_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[1] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G3xy_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[2] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G3xz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[3] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2x2y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[4] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2xyz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[5] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2x2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[6] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gx3y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[7] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[8] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[9] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gx3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[10] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G4y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[11] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G3yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[12] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2y2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[13] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gy3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[14] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G4z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[15] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G4x_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[16] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G3xy_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[17] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G3xz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[18] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2x2y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[19] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2xyz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[20] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2x2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[21] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gx3y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[22] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[23] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[24] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gx3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[25] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G4y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[26] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G3yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[27] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2y2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[28] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gy3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[29] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G4z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[30] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G4x_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[31] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G3xy_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[32] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G3xz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[33] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[34] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[35] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[36] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[37] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[38] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[39] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[40] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G4y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[41] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G3yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[42] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[43] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_a[44] += SQ_ERI_G_S_S_P_C1001003003_a_coefs*I_ERI_G4z_S_S_Pz_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_S_C3003_b
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_S_C3003_b_coefs = ic2*jc2*beta;
  SQ_ERI_G_S_S_S_C3003_b[0] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_G4x_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[1] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_G3xy_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[2] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_G3xz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[3] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_G2x2y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[4] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_G2xyz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[5] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_G2x2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[6] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_Gx3y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[7] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_Gx2yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[8] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_Gxy2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[9] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_Gx3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[10] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_G4y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[11] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_G3yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[12] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_G2y2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[13] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_Gy3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_b[14] += SQ_ERI_G_S_S_S_C3003_b_coefs*I_ERI_G4z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3003_b
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_S_S_C3003_b_coefs = ic2*jc2*beta;
  SQ_ERI_F_S_S_S_C3003_b[0] += SQ_ERI_F_S_S_S_C3003_b_coefs*I_ERI_F3x_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_b[1] += SQ_ERI_F_S_S_S_C3003_b_coefs*I_ERI_F2xy_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_b[2] += SQ_ERI_F_S_S_S_C3003_b_coefs*I_ERI_F2xz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_b[3] += SQ_ERI_F_S_S_S_C3003_b_coefs*I_ERI_Fx2y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_b[4] += SQ_ERI_F_S_S_S_C3003_b_coefs*I_ERI_Fxyz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_b[5] += SQ_ERI_F_S_S_S_C3003_b_coefs*I_ERI_Fx2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_b[6] += SQ_ERI_F_S_S_S_C3003_b_coefs*I_ERI_F3y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_b[7] += SQ_ERI_F_S_S_S_C3003_b_coefs*I_ERI_F2yz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_b[8] += SQ_ERI_F_S_S_S_C3003_b_coefs*I_ERI_Fy2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_b[9] += SQ_ERI_F_S_S_S_C3003_b_coefs*I_ERI_F3z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1003003_b
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1003003_b_coefs = ic2*jc2_1*beta;
  SQ_ERI_G_S_P_S_C1003003_b[0] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[1] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[2] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[3] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[4] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[5] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[6] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[7] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[8] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[9] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[10] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[11] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[12] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[13] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[14] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[15] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[16] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[17] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[18] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[19] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[20] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[21] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[22] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[23] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[24] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[25] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[26] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[27] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[28] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[29] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[30] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[31] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[32] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[33] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[34] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[35] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[36] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[37] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[38] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[39] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[40] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[41] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[42] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[43] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_b[44] += SQ_ERI_G_S_P_S_C1003003_b_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1003003_b
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1003003_b_coefs = ic2*jc2_1*beta;
  SQ_ERI_F_S_P_S_C1003003_b[0] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[1] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[2] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[3] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[4] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[5] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[6] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[7] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[8] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[9] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[10] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[11] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[12] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[13] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[14] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[15] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[16] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[17] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[18] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[19] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[20] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[21] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[22] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[23] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[24] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[25] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[26] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[27] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[28] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_b[29] += SQ_ERI_F_S_P_S_C1003003_b_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_P_C1000003003_b
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_P_C1000003003_b_coefs = ic2*jc2_2*beta;
  SQ_ERI_G_S_S_P_C1000003003_b[0] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G4x_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[1] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G3xy_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[2] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G3xz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[3] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2x2y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[4] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2xyz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[5] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2x2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[6] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gx3y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[7] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[8] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[9] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gx3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[10] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G4y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[11] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G3yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[12] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2y2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[13] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gy3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[14] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G4z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[15] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G4x_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[16] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G3xy_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[17] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G3xz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[18] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2x2y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[19] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2xyz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[20] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2x2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[21] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gx3y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[22] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[23] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[24] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gx3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[25] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G4y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[26] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G3yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[27] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2y2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[28] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gy3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[29] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G4z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[30] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G4x_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[31] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G3xy_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[32] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G3xz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[33] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[34] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[35] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[36] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[37] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[38] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[39] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[40] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G4y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[41] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G3yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[42] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[43] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_b[44] += SQ_ERI_G_S_S_P_C1000003003_b_coefs*I_ERI_G4z_S_S_Pz_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_S_C3003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_S_C3003_c_coefs = ic2*jc2*gamma;
  SQ_ERI_G_S_S_S_C3003_c[0] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_G4x_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[1] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_G3xy_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[2] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_G3xz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[3] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_G2x2y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[4] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_G2xyz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[5] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_G2x2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[6] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_Gx3y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[7] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_Gx2yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[8] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_Gxy2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[9] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_Gx3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[10] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_G4y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[11] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_G3yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[12] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_G2y2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[13] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_Gy3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_c[14] += SQ_ERI_G_S_S_S_C3003_c_coefs*I_ERI_G4z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_S_S_C3003_c_coefs = ic2*jc2*gamma;
  SQ_ERI_F_S_S_S_C3003_c[0] += SQ_ERI_F_S_S_S_C3003_c_coefs*I_ERI_F3x_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_c[1] += SQ_ERI_F_S_S_S_C3003_c_coefs*I_ERI_F2xy_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_c[2] += SQ_ERI_F_S_S_S_C3003_c_coefs*I_ERI_F2xz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_c[3] += SQ_ERI_F_S_S_S_C3003_c_coefs*I_ERI_Fx2y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_c[4] += SQ_ERI_F_S_S_S_C3003_c_coefs*I_ERI_Fxyz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_c[5] += SQ_ERI_F_S_S_S_C3003_c_coefs*I_ERI_Fx2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_c[6] += SQ_ERI_F_S_S_S_C3003_c_coefs*I_ERI_F3y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_c[7] += SQ_ERI_F_S_S_S_C3003_c_coefs*I_ERI_F2yz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_c[8] += SQ_ERI_F_S_S_S_C3003_c_coefs*I_ERI_Fy2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_c[9] += SQ_ERI_F_S_S_S_C3003_c_coefs*I_ERI_F3z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1003003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1003003_c_coefs = ic2*jc2_1*gamma;
  SQ_ERI_G_S_P_S_C1003003_c[0] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[1] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[2] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[3] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[4] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[5] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[6] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[7] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[8] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[9] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[10] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[11] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[12] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[13] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[14] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[15] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[16] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[17] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[18] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[19] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[20] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[21] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[22] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[23] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[24] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[25] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[26] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[27] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[28] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[29] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[30] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[31] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[32] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[33] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[34] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[35] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[36] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[37] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[38] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[39] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[40] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[41] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[42] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[43] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_c[44] += SQ_ERI_G_S_P_S_C1003003_c_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1003003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1003003_c_coefs = ic2*jc2_1*gamma;
  SQ_ERI_F_S_P_S_C1003003_c[0] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[1] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[2] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[3] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[4] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[5] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[6] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[7] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[8] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[9] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[10] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[11] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[12] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[13] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[14] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[15] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[16] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[17] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[18] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[19] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[20] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[21] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[22] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[23] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[24] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[25] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[26] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[27] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[28] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_c[29] += SQ_ERI_F_S_P_S_C1003003_c_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_P_C1000003003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_P_C1000003003_c_coefs = ic2*jc2_2*gamma;
  SQ_ERI_G_S_S_P_C1000003003_c[0] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G4x_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[1] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G3xy_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[2] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G3xz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[3] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2x2y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[4] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2xyz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[5] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2x2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[6] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gx3y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[7] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[8] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[9] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gx3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[10] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G4y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[11] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G3yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[12] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2y2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[13] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gy3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[14] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G4z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[15] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G4x_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[16] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G3xy_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[17] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G3xz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[18] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2x2y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[19] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2xyz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[20] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2x2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[21] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gx3y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[22] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[23] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[24] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gx3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[25] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G4y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[26] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G3yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[27] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2y2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[28] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gy3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[29] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G4z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[30] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G4x_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[31] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G3xy_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[32] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G3xz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[33] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[34] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[35] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[36] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[37] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[38] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[39] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[40] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G4y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[41] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G3yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[42] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[43] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_c[44] += SQ_ERI_G_S_S_P_C1000003003_c_coefs*I_ERI_G4z_S_S_Pz_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_S_C3003_ab
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_S_C3003_ab_coefs = ic2*jc2*alpha*beta;
  SQ_ERI_G_S_S_S_C3003_ab[0] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_G4x_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[1] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_G3xy_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[2] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_G3xz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[3] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_G2x2y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[4] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_G2xyz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[5] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_G2x2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[6] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_Gx3y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[7] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_Gx2yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[8] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_Gxy2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[9] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_Gx3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[10] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_G4y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[11] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_G3yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[12] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_G2y2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[13] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_Gy3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_ab[14] += SQ_ERI_G_S_S_S_C3003_ab_coefs*I_ERI_G4z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1003003_ab
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1003003_ab_coefs = ic2*jc2_1*alpha*beta;
  SQ_ERI_G_S_P_S_C1003003_ab[0] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[1] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[2] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[3] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[4] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[5] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[6] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[7] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[8] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[9] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[10] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[11] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[12] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[13] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[14] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[15] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[16] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[17] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[18] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[19] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[20] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[21] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[22] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[23] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[24] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[25] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[26] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[27] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[28] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[29] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[30] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[31] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[32] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[33] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[34] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[35] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[36] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[37] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[38] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[39] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[40] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[41] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[42] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[43] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_ab[44] += SQ_ERI_G_S_P_S_C1003003_ab_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_P_C1000003003_ab
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_P_C1000003003_ab_coefs = ic2*jc2_2*alpha*beta;
  SQ_ERI_G_S_S_P_C1000003003_ab[0] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G4x_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[1] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G3xy_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[2] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G3xz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[3] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2x2y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[4] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2xyz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[5] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2x2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[6] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gx3y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[7] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[8] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[9] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gx3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[10] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G4y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[11] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G3yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[12] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2y2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[13] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gy3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[14] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G4z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[15] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G4x_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[16] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G3xy_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[17] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G3xz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[18] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2x2y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[19] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2xyz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[20] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2x2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[21] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gx3y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[22] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[23] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[24] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gx3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[25] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G4y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[26] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G3yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[27] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2y2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[28] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gy3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[29] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G4z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[30] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G4x_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[31] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G3xy_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[32] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G3xz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[33] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[34] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[35] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[36] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[37] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[38] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[39] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[40] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G4y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[41] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G3yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[42] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[43] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_ab[44] += SQ_ERI_G_S_S_P_C1000003003_ab_coefs*I_ERI_G4z_S_S_Pz_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C3003_bc
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C3003_bc_coefs = ic2*jc2*beta*gamma;
  SQ_ERI_G_S_P_S_C3003_bc[0] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[1] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[2] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[3] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[4] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[5] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[6] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[7] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[8] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[9] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[10] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[11] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[12] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[13] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[14] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[15] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[16] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[17] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[18] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[19] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[20] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[21] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[22] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[23] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[24] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[25] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[26] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[27] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[28] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[29] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[30] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[31] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[32] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[33] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[34] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[35] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[36] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[37] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[38] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[39] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[40] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[41] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[42] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[43] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C3003_bc[44] += SQ_ERI_G_S_P_S_C3003_bc_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C3003_bc
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C3003_bc_coefs = ic2*jc2*beta*gamma;
  SQ_ERI_F_S_P_S_C3003_bc[0] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[1] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[2] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[3] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[4] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[5] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[6] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[7] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[8] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[9] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[10] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[11] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[12] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[13] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[14] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[15] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[16] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[17] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[18] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[19] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[20] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[21] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[22] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[23] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[24] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[25] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[26] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[27] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[28] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C3003_bc[29] += SQ_ERI_F_S_P_S_C3003_bc_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_S_C1003003_b
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_S_C1003003_b_coefs = ic2*jc2_1*beta;
  SQ_ERI_G_S_S_S_C1003003_b[0] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_G4x_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[1] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_G3xy_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[2] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_G3xz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[3] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_G2x2y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[4] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_G2xyz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[5] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_G2x2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[6] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_Gx3y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[7] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_Gx2yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[8] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_Gxy2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[9] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_Gx3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[10] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_G4y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[11] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_G3yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[12] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_G2y2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[13] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_Gy3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C1003003_b[14] += SQ_ERI_G_S_S_S_C1003003_b_coefs*I_ERI_G4z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C1003003_b
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_S_S_C1003003_b_coefs = ic2*jc2_1*beta;
  SQ_ERI_F_S_S_S_C1003003_b[0] += SQ_ERI_F_S_S_S_C1003003_b_coefs*I_ERI_F3x_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003_b[1] += SQ_ERI_F_S_S_S_C1003003_b_coefs*I_ERI_F2xy_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003_b[2] += SQ_ERI_F_S_S_S_C1003003_b_coefs*I_ERI_F2xz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003_b[3] += SQ_ERI_F_S_S_S_C1003003_b_coefs*I_ERI_Fx2y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003_b[4] += SQ_ERI_F_S_S_S_C1003003_b_coefs*I_ERI_Fxyz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003_b[5] += SQ_ERI_F_S_S_S_C1003003_b_coefs*I_ERI_Fx2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003_b[6] += SQ_ERI_F_S_S_S_C1003003_b_coefs*I_ERI_F3y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003_b[7] += SQ_ERI_F_S_S_S_C1003003_b_coefs*I_ERI_F2yz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003_b[8] += SQ_ERI_F_S_S_S_C1003003_b_coefs*I_ERI_Fy2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C1003003_b[9] += SQ_ERI_F_S_S_S_C1003003_b_coefs*I_ERI_F3z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_P_C1001003003_b
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_P_C1001003003_b_coefs = ic2*jc2_3*beta;
  SQ_ERI_G_S_S_P_C1001003003_b[0] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G4x_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[1] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G3xy_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[2] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G3xz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[3] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2x2y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[4] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2xyz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[5] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2x2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[6] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gx3y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[7] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[8] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[9] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gx3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[10] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G4y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[11] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G3yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[12] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2y2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[13] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gy3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[14] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G4z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[15] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G4x_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[16] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G3xy_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[17] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G3xz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[18] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2x2y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[19] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2xyz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[20] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2x2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[21] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gx3y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[22] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[23] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[24] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gx3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[25] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G4y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[26] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G3yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[27] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2y2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[28] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gy3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[29] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G4z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[30] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G4x_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[31] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G3xy_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[32] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G3xz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[33] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[34] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[35] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[36] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[37] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[38] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[39] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[40] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G4y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[41] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G3yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[42] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[43] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1001003003_b[44] += SQ_ERI_G_S_S_P_C1001003003_b_coefs*I_ERI_G4z_S_S_Pz_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_S_C3003_bb
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_S_C3003_bb_coefs = ic2*jc2*beta*beta;
  SQ_ERI_G_S_S_S_C3003_bb[0] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_G4x_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[1] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_G3xy_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[2] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_G3xz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[3] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_G2x2y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[4] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_G2xyz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[5] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_G2x2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[6] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_Gx3y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[7] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_Gx2yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[8] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_Gxy2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[9] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_Gx3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[10] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_G4y_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[11] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_G3yz_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[12] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_G2y2z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[13] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_Gy3z_S_S_S_vrr;
  SQ_ERI_G_S_S_S_C3003_bb[14] += SQ_ERI_G_S_S_S_C3003_bb_coefs*I_ERI_G4z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_S_S_C3003_bb
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_S_S_C3003_bb_coefs = ic2*jc2*beta*beta;
  SQ_ERI_F_S_S_S_C3003_bb[0] += SQ_ERI_F_S_S_S_C3003_bb_coefs*I_ERI_F3x_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_bb[1] += SQ_ERI_F_S_S_S_C3003_bb_coefs*I_ERI_F2xy_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_bb[2] += SQ_ERI_F_S_S_S_C3003_bb_coefs*I_ERI_F2xz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_bb[3] += SQ_ERI_F_S_S_S_C3003_bb_coefs*I_ERI_Fx2y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_bb[4] += SQ_ERI_F_S_S_S_C3003_bb_coefs*I_ERI_Fxyz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_bb[5] += SQ_ERI_F_S_S_S_C3003_bb_coefs*I_ERI_Fx2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_bb[6] += SQ_ERI_F_S_S_S_C3003_bb_coefs*I_ERI_F3y_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_bb[7] += SQ_ERI_F_S_S_S_C3003_bb_coefs*I_ERI_F2yz_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_bb[8] += SQ_ERI_F_S_S_S_C3003_bb_coefs*I_ERI_Fy2z_S_S_S_vrr;
  SQ_ERI_F_S_S_S_C3003_bb[9] += SQ_ERI_F_S_S_S_C3003_bb_coefs*I_ERI_F3z_S_S_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1003003_bb
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1003003_bb_coefs = ic2*jc2_1*beta*beta;
  SQ_ERI_G_S_P_S_C1003003_bb[0] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[1] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[2] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[3] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[4] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[5] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[6] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[7] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[8] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[9] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[10] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[11] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[12] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[13] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[14] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[15] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[16] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[17] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[18] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[19] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[20] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[21] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[22] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[23] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[24] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[25] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[26] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[27] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[28] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[29] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[30] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[31] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[32] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[33] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[34] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[35] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[36] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[37] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[38] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[39] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[40] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[41] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[42] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[43] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1003003_bb[44] += SQ_ERI_G_S_P_S_C1003003_bb_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1003003_bb
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1003003_bb_coefs = ic2*jc2_1*beta*beta;
  SQ_ERI_F_S_P_S_C1003003_bb[0] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[1] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[2] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[3] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[4] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[5] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[6] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[7] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[8] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[9] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[10] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[11] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[12] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[13] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[14] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[15] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[16] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[17] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[18] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[19] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[20] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[21] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[22] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[23] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[24] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[25] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[26] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[27] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[28] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1003003_bb[29] += SQ_ERI_F_S_P_S_C1003003_bb_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_S_P_C1000003003_bb
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_S_P_C1000003003_bb_coefs = ic2*jc2_2*beta*beta;
  SQ_ERI_G_S_S_P_C1000003003_bb[0] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G4x_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[1] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G3xy_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[2] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G3xz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[3] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2x2y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[4] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2xyz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[5] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2x2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[6] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gx3y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[7] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gx2yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[8] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gxy2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[9] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gx3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[10] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G4y_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[11] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G3yz_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[12] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2y2z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[13] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gy3z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[14] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G4z_S_S_Px_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[15] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G4x_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[16] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G3xy_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[17] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G3xz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[18] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2x2y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[19] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2xyz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[20] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2x2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[21] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gx3y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[22] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gx2yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[23] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gxy2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[24] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gx3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[25] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G4y_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[26] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G3yz_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[27] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2y2z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[28] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gy3z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[29] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G4z_S_S_Py_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[30] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G4x_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[31] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G3xy_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[32] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G3xz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[33] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2x2y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[34] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2xyz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[35] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2x2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[36] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gx3y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[37] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gx2yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[38] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gxy2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[39] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gx3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[40] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G4y_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[41] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G3yz_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[42] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G2y2z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[43] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_Gy3z_S_S_Pz_vrr;
  SQ_ERI_G_S_S_P_C1000003003_bb[44] += SQ_ERI_G_S_S_P_C1000003003_bb_coefs*I_ERI_G4z_S_S_Pz_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1001003003
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1001003003_coefs = ic2*jc2_3;
  SQ_ERI_G_S_P_S_C1001003003[0] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[1] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[2] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[3] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[4] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[5] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[6] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[7] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[8] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[9] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[10] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[11] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[12] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[13] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[14] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[15] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[16] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[17] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[18] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[19] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[20] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[21] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[22] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[23] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[24] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[25] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[26] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[27] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[28] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[29] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[30] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[31] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[32] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[33] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[34] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[35] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[36] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[37] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[38] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[39] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[40] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[41] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[42] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[43] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003[44] += SQ_ERI_G_S_P_S_C1001003003_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1001003003
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1001003003_coefs = ic2*jc2_3;
  SQ_ERI_F_S_P_S_C1001003003[0] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[1] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[2] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[3] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[4] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[5] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[6] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[7] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[8] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[9] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[10] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[11] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[12] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[13] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[14] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[15] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[16] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[17] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[18] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[19] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[20] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[21] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[22] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[23] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[24] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[25] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[26] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[27] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[28] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003[29] += SQ_ERI_F_S_P_S_C1001003003_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1001003003_a
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1001003003_a_coefs = ic2*jc2_3*alpha;
  SQ_ERI_G_S_P_S_C1001003003_a[0] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[1] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[2] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[3] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[4] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[5] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[6] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[7] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[8] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[9] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[10] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[11] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[12] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[13] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[14] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[15] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[16] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[17] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[18] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[19] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[20] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[21] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[22] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[23] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[24] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[25] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[26] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[27] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[28] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[29] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[30] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[31] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[32] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[33] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[34] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[35] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[36] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[37] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[38] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[39] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[40] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[41] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[42] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[43] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_a[44] += SQ_ERI_G_S_P_S_C1001003003_a_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1000003003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1000003003_c_coefs = ic2*jc2_2*gamma;
  SQ_ERI_G_S_P_S_C1000003003_c[0] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[1] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[2] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[3] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[4] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[5] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[6] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[7] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[8] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[9] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[10] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[11] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[12] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[13] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[14] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[15] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[16] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[17] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[18] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[19] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[20] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[21] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[22] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[23] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[24] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[25] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[26] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[27] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[28] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[29] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[30] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[31] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[32] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[33] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[34] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[35] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[36] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[37] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[38] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[39] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[40] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[41] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[42] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[43] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_c[44] += SQ_ERI_G_S_P_S_C1000003003_c_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1000003003_c_coefs = ic2*jc2_2*gamma;
  SQ_ERI_F_S_P_S_C1000003003_c[0] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[1] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[2] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[3] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[4] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[5] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[6] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[7] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[8] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[9] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[10] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[11] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[12] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[13] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[14] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[15] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[16] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[17] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[18] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[19] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[20] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[21] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[22] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[23] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[24] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[25] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[26] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[27] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[28] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_c[29] += SQ_ERI_F_S_P_S_C1000003003_c_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1001003003_a
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1001003003_a_coefs = ic2*jc2_3*alpha;
  SQ_ERI_F_S_P_S_C1001003003_a[0] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[1] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[2] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[3] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[4] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[5] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[6] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[7] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[8] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[9] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[10] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[11] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[12] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[13] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[14] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[15] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[16] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[17] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[18] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[19] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[20] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[21] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[22] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[23] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[24] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[25] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[26] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[27] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[28] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_a[29] += SQ_ERI_F_S_P_S_C1001003003_a_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1000003003_ac
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1000003003_ac_coefs = ic2*jc2_2*alpha*gamma;
  SQ_ERI_G_S_P_S_C1000003003_ac[0] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[1] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[2] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[3] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[4] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[5] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[6] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[7] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[8] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[9] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[10] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[11] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[12] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[13] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[14] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[15] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[16] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[17] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[18] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[19] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[20] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[21] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[22] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[23] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[24] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[25] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[26] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[27] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[28] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[29] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[30] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[31] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[32] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[33] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[34] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[35] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[36] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[37] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[38] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[39] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[40] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[41] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[42] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[43] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_ac[44] += SQ_ERI_G_S_P_S_C1000003003_ac_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1001003003_b
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1001003003_b_coefs = ic2*jc2_3*beta;
  SQ_ERI_G_S_P_S_C1001003003_b[0] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[1] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[2] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[3] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[4] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[5] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[6] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[7] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[8] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[9] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[10] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[11] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[12] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[13] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[14] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[15] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[16] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[17] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[18] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[19] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[20] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[21] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[22] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[23] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[24] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[25] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[26] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[27] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[28] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[29] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[30] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[31] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[32] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[33] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[34] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[35] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[36] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[37] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[38] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[39] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[40] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[41] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[42] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[43] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_b[44] += SQ_ERI_G_S_P_S_C1001003003_b_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1001003003_b
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1001003003_b_coefs = ic2*jc2_3*beta;
  SQ_ERI_F_S_P_S_C1001003003_b[0] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[1] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[2] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[3] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[4] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[5] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[6] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[7] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[8] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[9] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[10] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[11] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[12] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[13] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[14] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[15] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[16] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[17] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[18] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[19] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[20] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[21] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[22] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[23] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[24] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[25] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[26] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[27] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[28] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_b[29] += SQ_ERI_F_S_P_S_C1001003003_b_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1001003003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1001003003_c_coefs = ic2*jc2_3*gamma;
  SQ_ERI_G_S_P_S_C1001003003_c[0] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[1] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[2] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[3] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[4] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[5] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[6] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[7] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[8] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[9] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[10] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[11] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[12] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[13] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[14] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[15] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[16] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[17] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[18] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[19] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[20] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[21] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[22] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[23] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[24] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[25] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[26] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[27] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[28] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[29] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[30] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[31] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[32] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[33] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[34] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[35] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[36] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[37] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[38] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[39] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[40] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[41] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[42] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[43] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_c[44] += SQ_ERI_G_S_P_S_C1001003003_c_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1001003003_c
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1001003003_c_coefs = ic2*jc2_3*gamma;
  SQ_ERI_F_S_P_S_C1001003003_c[0] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[1] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[2] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[3] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[4] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[5] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[6] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[7] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[8] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[9] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[10] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[11] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[12] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[13] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[14] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[15] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[16] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[17] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[18] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[19] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[20] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[21] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[22] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[23] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[24] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[25] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[26] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[27] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[28] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_c[29] += SQ_ERI_F_S_P_S_C1001003003_c_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1001003003_ab
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1001003003_ab_coefs = ic2*jc2_3*alpha*beta;
  SQ_ERI_G_S_P_S_C1001003003_ab[0] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[1] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[2] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[3] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[4] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[5] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[6] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[7] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[8] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[9] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[10] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[11] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[12] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[13] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[14] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[15] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[16] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[17] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[18] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[19] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[20] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[21] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[22] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[23] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[24] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[25] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[26] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[27] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[28] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[29] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[30] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[31] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[32] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[33] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[34] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[35] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[36] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[37] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[38] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[39] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[40] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[41] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[42] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[43] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_ab[44] += SQ_ERI_G_S_P_S_C1001003003_ab_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1000003003_bc
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1000003003_bc_coefs = ic2*jc2_2*beta*gamma;
  SQ_ERI_G_S_P_S_C1000003003_bc[0] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[1] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[2] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[3] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[4] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[5] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[6] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[7] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[8] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[9] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[10] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[11] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[12] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[13] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[14] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[15] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[16] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[17] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[18] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[19] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[20] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[21] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[22] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[23] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[24] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[25] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[26] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[27] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[28] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[29] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[30] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[31] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[32] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[33] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[34] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[35] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[36] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[37] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[38] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[39] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[40] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[41] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[42] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[43] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1000003003_bc[44] += SQ_ERI_G_S_P_S_C1000003003_bc_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1000003003_bc
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1000003003_bc_coefs = ic2*jc2_2*beta*gamma;
  SQ_ERI_F_S_P_S_C1000003003_bc[0] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[1] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[2] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[3] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[4] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[5] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[6] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[7] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[8] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[9] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[10] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[11] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[12] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[13] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[14] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[15] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[16] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[17] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[18] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[19] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[20] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[21] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[22] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[23] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[24] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[25] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[26] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[27] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[28] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1000003003_bc[29] += SQ_ERI_F_S_P_S_C1000003003_bc_coefs*I_ERI_F3z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_S_C1001003003_bb
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_P_S_C1001003003_bb_coefs = ic2*jc2_3*beta*beta;
  SQ_ERI_G_S_P_S_C1001003003_bb[0] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G4x_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[1] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G3xy_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[2] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G3xz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[3] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2x2y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[4] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2xyz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[5] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2x2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[6] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gx3y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[7] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[8] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[9] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gx3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[10] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G4y_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[11] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G3yz_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[12] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2y2z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[13] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gy3z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[14] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G4z_S_Px_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[15] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G4x_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[16] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G3xy_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[17] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G3xz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[18] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2x2y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[19] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2xyz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[20] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2x2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[21] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gx3y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[22] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[23] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[24] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gx3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[25] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G4y_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[26] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G3yz_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[27] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2y2z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[28] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gy3z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[29] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G4z_S_Py_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[30] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G4x_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[31] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G3xy_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[32] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G3xz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[33] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[34] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[35] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[36] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[37] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[38] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[39] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[40] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G4y_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[41] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G3yz_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[42] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[43] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
  SQ_ERI_G_S_P_S_C1001003003_bb[44] += SQ_ERI_G_S_P_S_C1001003003_bb_coefs*I_ERI_G4z_S_Pz_S_vrr;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_P_S_C1001003003_bb
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_F_S_P_S_C1001003003_bb_coefs = ic2*jc2_3*beta*beta;
  SQ_ERI_F_S_P_S_C1001003003_bb[0] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F3x_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[1] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F2xy_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[2] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F2xz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[3] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fx2y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[4] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fxyz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[5] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fx2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[6] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F3y_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[7] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F2yz_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[8] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fy2z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[9] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F3z_S_Px_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[10] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F3x_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[11] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F2xy_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[12] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F2xz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[13] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fx2y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[14] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fxyz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[15] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fx2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[16] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F3y_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[17] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F2yz_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[18] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fy2z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[19] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F3z_S_Py_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[20] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F3x_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[21] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F2xy_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[22] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F2xz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[23] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[24] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[25] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[26] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F3y_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[27] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F2yz_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[28] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
  SQ_ERI_F_S_P_S_C1001003003_bb[29] += SQ_ERI_F_S_P_S_C1001003003_bb_coefs*I_ERI_F3z_S_Pz_S_vrr;
}
