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
// BRA1 as redundant position, total RHS integrals evaluated as: 614854
// BRA2 as redundant position, total RHS integrals evaluated as: 587361
// KET1 as redundant position, total RHS integrals evaluated as: 614845
// KET2 as redundant position, total RHS integrals evaluated as: 587361
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

void hgp_os_eri_f_sp_f_sp_d2_vrrcont_4(const Double& ic2, const Double& ic2_1, const Double& jc2, const Double& jc2_1, const Double& alpha, const Double& beta, const Double& gamma, const Double& delta, const Double* SQ_ERI_G_S_I_S_vrr_array, Double* SQ_ERI_G_S_I_S_C1003001003_cc  ) 
{

  /************************************************************
   * transform the array form of integral into variables: SQ_ERI_G_S_I_S
   ************************************************************/
  Double I_ERI_G4x_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[0];
  Double I_ERI_G3xy_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[1];
  Double I_ERI_G3xz_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[2];
  Double I_ERI_G2x2y_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[3];
  Double I_ERI_G2xyz_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[4];
  Double I_ERI_G2x2z_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[5];
  Double I_ERI_Gx3y_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[6];
  Double I_ERI_Gx2yz_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[7];
  Double I_ERI_Gxy2z_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[8];
  Double I_ERI_Gx3z_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[9];
  Double I_ERI_G4y_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[10];
  Double I_ERI_G3yz_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[11];
  Double I_ERI_G2y2z_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[12];
  Double I_ERI_Gy3z_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[13];
  Double I_ERI_G4z_S_I6x_S_vrr = SQ_ERI_G_S_I_S_vrr_array[14];
  Double I_ERI_G4x_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[15];
  Double I_ERI_G3xy_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[16];
  Double I_ERI_G3xz_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[17];
  Double I_ERI_G2x2y_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[18];
  Double I_ERI_G2xyz_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[19];
  Double I_ERI_G2x2z_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[20];
  Double I_ERI_Gx3y_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[21];
  Double I_ERI_Gx2yz_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[22];
  Double I_ERI_Gxy2z_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[23];
  Double I_ERI_Gx3z_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[24];
  Double I_ERI_G4y_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[25];
  Double I_ERI_G3yz_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[26];
  Double I_ERI_G2y2z_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[27];
  Double I_ERI_Gy3z_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[28];
  Double I_ERI_G4z_S_I5xy_S_vrr = SQ_ERI_G_S_I_S_vrr_array[29];
  Double I_ERI_G4x_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[30];
  Double I_ERI_G3xy_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[31];
  Double I_ERI_G3xz_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[32];
  Double I_ERI_G2x2y_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[33];
  Double I_ERI_G2xyz_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[34];
  Double I_ERI_G2x2z_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[35];
  Double I_ERI_Gx3y_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[36];
  Double I_ERI_Gx2yz_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[37];
  Double I_ERI_Gxy2z_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[38];
  Double I_ERI_Gx3z_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[39];
  Double I_ERI_G4y_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[40];
  Double I_ERI_G3yz_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[41];
  Double I_ERI_G2y2z_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[42];
  Double I_ERI_Gy3z_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[43];
  Double I_ERI_G4z_S_I5xz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[44];
  Double I_ERI_G4x_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[45];
  Double I_ERI_G3xy_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[46];
  Double I_ERI_G3xz_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[47];
  Double I_ERI_G2x2y_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[48];
  Double I_ERI_G2xyz_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[49];
  Double I_ERI_G2x2z_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[50];
  Double I_ERI_Gx3y_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[51];
  Double I_ERI_Gx2yz_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[52];
  Double I_ERI_Gxy2z_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[53];
  Double I_ERI_Gx3z_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[54];
  Double I_ERI_G4y_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[55];
  Double I_ERI_G3yz_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[56];
  Double I_ERI_G2y2z_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[57];
  Double I_ERI_Gy3z_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[58];
  Double I_ERI_G4z_S_I4x2y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[59];
  Double I_ERI_G4x_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[60];
  Double I_ERI_G3xy_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[61];
  Double I_ERI_G3xz_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[62];
  Double I_ERI_G2x2y_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[63];
  Double I_ERI_G2xyz_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[64];
  Double I_ERI_G2x2z_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[65];
  Double I_ERI_Gx3y_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[66];
  Double I_ERI_Gx2yz_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[67];
  Double I_ERI_Gxy2z_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[68];
  Double I_ERI_Gx3z_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[69];
  Double I_ERI_G4y_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[70];
  Double I_ERI_G3yz_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[71];
  Double I_ERI_G2y2z_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[72];
  Double I_ERI_Gy3z_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[73];
  Double I_ERI_G4z_S_I4xyz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[74];
  Double I_ERI_G4x_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[75];
  Double I_ERI_G3xy_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[76];
  Double I_ERI_G3xz_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[77];
  Double I_ERI_G2x2y_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[78];
  Double I_ERI_G2xyz_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[79];
  Double I_ERI_G2x2z_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[80];
  Double I_ERI_Gx3y_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[81];
  Double I_ERI_Gx2yz_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[82];
  Double I_ERI_Gxy2z_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[83];
  Double I_ERI_Gx3z_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[84];
  Double I_ERI_G4y_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[85];
  Double I_ERI_G3yz_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[86];
  Double I_ERI_G2y2z_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[87];
  Double I_ERI_Gy3z_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[88];
  Double I_ERI_G4z_S_I4x2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[89];
  Double I_ERI_G4x_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[90];
  Double I_ERI_G3xy_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[91];
  Double I_ERI_G3xz_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[92];
  Double I_ERI_G2x2y_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[93];
  Double I_ERI_G2xyz_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[94];
  Double I_ERI_G2x2z_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[95];
  Double I_ERI_Gx3y_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[96];
  Double I_ERI_Gx2yz_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[97];
  Double I_ERI_Gxy2z_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[98];
  Double I_ERI_Gx3z_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[99];
  Double I_ERI_G4y_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[100];
  Double I_ERI_G3yz_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[101];
  Double I_ERI_G2y2z_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[102];
  Double I_ERI_Gy3z_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[103];
  Double I_ERI_G4z_S_I3x3y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[104];
  Double I_ERI_G4x_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[105];
  Double I_ERI_G3xy_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[106];
  Double I_ERI_G3xz_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[107];
  Double I_ERI_G2x2y_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[108];
  Double I_ERI_G2xyz_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[109];
  Double I_ERI_G2x2z_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[110];
  Double I_ERI_Gx3y_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[111];
  Double I_ERI_Gx2yz_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[112];
  Double I_ERI_Gxy2z_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[113];
  Double I_ERI_Gx3z_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[114];
  Double I_ERI_G4y_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[115];
  Double I_ERI_G3yz_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[116];
  Double I_ERI_G2y2z_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[117];
  Double I_ERI_Gy3z_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[118];
  Double I_ERI_G4z_S_I3x2yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[119];
  Double I_ERI_G4x_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[120];
  Double I_ERI_G3xy_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[121];
  Double I_ERI_G3xz_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[122];
  Double I_ERI_G2x2y_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[123];
  Double I_ERI_G2xyz_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[124];
  Double I_ERI_G2x2z_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[125];
  Double I_ERI_Gx3y_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[126];
  Double I_ERI_Gx2yz_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[127];
  Double I_ERI_Gxy2z_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[128];
  Double I_ERI_Gx3z_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[129];
  Double I_ERI_G4y_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[130];
  Double I_ERI_G3yz_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[131];
  Double I_ERI_G2y2z_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[132];
  Double I_ERI_Gy3z_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[133];
  Double I_ERI_G4z_S_I3xy2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[134];
  Double I_ERI_G4x_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[135];
  Double I_ERI_G3xy_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[136];
  Double I_ERI_G3xz_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[137];
  Double I_ERI_G2x2y_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[138];
  Double I_ERI_G2xyz_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[139];
  Double I_ERI_G2x2z_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[140];
  Double I_ERI_Gx3y_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[141];
  Double I_ERI_Gx2yz_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[142];
  Double I_ERI_Gxy2z_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[143];
  Double I_ERI_Gx3z_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[144];
  Double I_ERI_G4y_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[145];
  Double I_ERI_G3yz_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[146];
  Double I_ERI_G2y2z_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[147];
  Double I_ERI_Gy3z_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[148];
  Double I_ERI_G4z_S_I3x3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[149];
  Double I_ERI_G4x_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[150];
  Double I_ERI_G3xy_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[151];
  Double I_ERI_G3xz_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[152];
  Double I_ERI_G2x2y_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[153];
  Double I_ERI_G2xyz_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[154];
  Double I_ERI_G2x2z_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[155];
  Double I_ERI_Gx3y_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[156];
  Double I_ERI_Gx2yz_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[157];
  Double I_ERI_Gxy2z_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[158];
  Double I_ERI_Gx3z_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[159];
  Double I_ERI_G4y_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[160];
  Double I_ERI_G3yz_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[161];
  Double I_ERI_G2y2z_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[162];
  Double I_ERI_Gy3z_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[163];
  Double I_ERI_G4z_S_I2x4y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[164];
  Double I_ERI_G4x_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[165];
  Double I_ERI_G3xy_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[166];
  Double I_ERI_G3xz_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[167];
  Double I_ERI_G2x2y_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[168];
  Double I_ERI_G2xyz_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[169];
  Double I_ERI_G2x2z_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[170];
  Double I_ERI_Gx3y_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[171];
  Double I_ERI_Gx2yz_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[172];
  Double I_ERI_Gxy2z_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[173];
  Double I_ERI_Gx3z_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[174];
  Double I_ERI_G4y_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[175];
  Double I_ERI_G3yz_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[176];
  Double I_ERI_G2y2z_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[177];
  Double I_ERI_Gy3z_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[178];
  Double I_ERI_G4z_S_I2x3yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[179];
  Double I_ERI_G4x_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[180];
  Double I_ERI_G3xy_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[181];
  Double I_ERI_G3xz_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[182];
  Double I_ERI_G2x2y_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[183];
  Double I_ERI_G2xyz_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[184];
  Double I_ERI_G2x2z_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[185];
  Double I_ERI_Gx3y_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[186];
  Double I_ERI_Gx2yz_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[187];
  Double I_ERI_Gxy2z_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[188];
  Double I_ERI_Gx3z_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[189];
  Double I_ERI_G4y_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[190];
  Double I_ERI_G3yz_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[191];
  Double I_ERI_G2y2z_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[192];
  Double I_ERI_Gy3z_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[193];
  Double I_ERI_G4z_S_I2x2y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[194];
  Double I_ERI_G4x_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[195];
  Double I_ERI_G3xy_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[196];
  Double I_ERI_G3xz_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[197];
  Double I_ERI_G2x2y_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[198];
  Double I_ERI_G2xyz_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[199];
  Double I_ERI_G2x2z_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[200];
  Double I_ERI_Gx3y_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[201];
  Double I_ERI_Gx2yz_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[202];
  Double I_ERI_Gxy2z_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[203];
  Double I_ERI_Gx3z_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[204];
  Double I_ERI_G4y_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[205];
  Double I_ERI_G3yz_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[206];
  Double I_ERI_G2y2z_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[207];
  Double I_ERI_Gy3z_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[208];
  Double I_ERI_G4z_S_I2xy3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[209];
  Double I_ERI_G4x_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[210];
  Double I_ERI_G3xy_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[211];
  Double I_ERI_G3xz_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[212];
  Double I_ERI_G2x2y_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[213];
  Double I_ERI_G2xyz_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[214];
  Double I_ERI_G2x2z_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[215];
  Double I_ERI_Gx3y_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[216];
  Double I_ERI_Gx2yz_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[217];
  Double I_ERI_Gxy2z_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[218];
  Double I_ERI_Gx3z_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[219];
  Double I_ERI_G4y_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[220];
  Double I_ERI_G3yz_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[221];
  Double I_ERI_G2y2z_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[222];
  Double I_ERI_Gy3z_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[223];
  Double I_ERI_G4z_S_I2x4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[224];
  Double I_ERI_G4x_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[225];
  Double I_ERI_G3xy_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[226];
  Double I_ERI_G3xz_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[227];
  Double I_ERI_G2x2y_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[228];
  Double I_ERI_G2xyz_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[229];
  Double I_ERI_G2x2z_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[230];
  Double I_ERI_Gx3y_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[231];
  Double I_ERI_Gx2yz_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[232];
  Double I_ERI_Gxy2z_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[233];
  Double I_ERI_Gx3z_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[234];
  Double I_ERI_G4y_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[235];
  Double I_ERI_G3yz_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[236];
  Double I_ERI_G2y2z_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[237];
  Double I_ERI_Gy3z_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[238];
  Double I_ERI_G4z_S_Ix5y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[239];
  Double I_ERI_G4x_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[240];
  Double I_ERI_G3xy_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[241];
  Double I_ERI_G3xz_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[242];
  Double I_ERI_G2x2y_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[243];
  Double I_ERI_G2xyz_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[244];
  Double I_ERI_G2x2z_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[245];
  Double I_ERI_Gx3y_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[246];
  Double I_ERI_Gx2yz_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[247];
  Double I_ERI_Gxy2z_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[248];
  Double I_ERI_Gx3z_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[249];
  Double I_ERI_G4y_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[250];
  Double I_ERI_G3yz_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[251];
  Double I_ERI_G2y2z_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[252];
  Double I_ERI_Gy3z_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[253];
  Double I_ERI_G4z_S_Ix4yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[254];
  Double I_ERI_G4x_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[255];
  Double I_ERI_G3xy_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[256];
  Double I_ERI_G3xz_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[257];
  Double I_ERI_G2x2y_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[258];
  Double I_ERI_G2xyz_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[259];
  Double I_ERI_G2x2z_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[260];
  Double I_ERI_Gx3y_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[261];
  Double I_ERI_Gx2yz_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[262];
  Double I_ERI_Gxy2z_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[263];
  Double I_ERI_Gx3z_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[264];
  Double I_ERI_G4y_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[265];
  Double I_ERI_G3yz_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[266];
  Double I_ERI_G2y2z_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[267];
  Double I_ERI_Gy3z_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[268];
  Double I_ERI_G4z_S_Ix3y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[269];
  Double I_ERI_G4x_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[270];
  Double I_ERI_G3xy_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[271];
  Double I_ERI_G3xz_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[272];
  Double I_ERI_G2x2y_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[273];
  Double I_ERI_G2xyz_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[274];
  Double I_ERI_G2x2z_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[275];
  Double I_ERI_Gx3y_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[276];
  Double I_ERI_Gx2yz_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[277];
  Double I_ERI_Gxy2z_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[278];
  Double I_ERI_Gx3z_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[279];
  Double I_ERI_G4y_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[280];
  Double I_ERI_G3yz_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[281];
  Double I_ERI_G2y2z_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[282];
  Double I_ERI_Gy3z_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[283];
  Double I_ERI_G4z_S_Ix2y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[284];
  Double I_ERI_G4x_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[285];
  Double I_ERI_G3xy_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[286];
  Double I_ERI_G3xz_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[287];
  Double I_ERI_G2x2y_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[288];
  Double I_ERI_G2xyz_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[289];
  Double I_ERI_G2x2z_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[290];
  Double I_ERI_Gx3y_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[291];
  Double I_ERI_Gx2yz_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[292];
  Double I_ERI_Gxy2z_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[293];
  Double I_ERI_Gx3z_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[294];
  Double I_ERI_G4y_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[295];
  Double I_ERI_G3yz_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[296];
  Double I_ERI_G2y2z_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[297];
  Double I_ERI_Gy3z_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[298];
  Double I_ERI_G4z_S_Ixy4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[299];
  Double I_ERI_G4x_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[300];
  Double I_ERI_G3xy_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[301];
  Double I_ERI_G3xz_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[302];
  Double I_ERI_G2x2y_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[303];
  Double I_ERI_G2xyz_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[304];
  Double I_ERI_G2x2z_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[305];
  Double I_ERI_Gx3y_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[306];
  Double I_ERI_Gx2yz_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[307];
  Double I_ERI_Gxy2z_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[308];
  Double I_ERI_Gx3z_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[309];
  Double I_ERI_G4y_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[310];
  Double I_ERI_G3yz_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[311];
  Double I_ERI_G2y2z_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[312];
  Double I_ERI_Gy3z_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[313];
  Double I_ERI_G4z_S_Ix5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[314];
  Double I_ERI_G4x_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[315];
  Double I_ERI_G3xy_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[316];
  Double I_ERI_G3xz_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[317];
  Double I_ERI_G2x2y_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[318];
  Double I_ERI_G2xyz_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[319];
  Double I_ERI_G2x2z_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[320];
  Double I_ERI_Gx3y_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[321];
  Double I_ERI_Gx2yz_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[322];
  Double I_ERI_Gxy2z_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[323];
  Double I_ERI_Gx3z_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[324];
  Double I_ERI_G4y_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[325];
  Double I_ERI_G3yz_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[326];
  Double I_ERI_G2y2z_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[327];
  Double I_ERI_Gy3z_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[328];
  Double I_ERI_G4z_S_I6y_S_vrr = SQ_ERI_G_S_I_S_vrr_array[329];
  Double I_ERI_G4x_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[330];
  Double I_ERI_G3xy_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[331];
  Double I_ERI_G3xz_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[332];
  Double I_ERI_G2x2y_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[333];
  Double I_ERI_G2xyz_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[334];
  Double I_ERI_G2x2z_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[335];
  Double I_ERI_Gx3y_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[336];
  Double I_ERI_Gx2yz_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[337];
  Double I_ERI_Gxy2z_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[338];
  Double I_ERI_Gx3z_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[339];
  Double I_ERI_G4y_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[340];
  Double I_ERI_G3yz_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[341];
  Double I_ERI_G2y2z_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[342];
  Double I_ERI_Gy3z_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[343];
  Double I_ERI_G4z_S_I5yz_S_vrr = SQ_ERI_G_S_I_S_vrr_array[344];
  Double I_ERI_G4x_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[345];
  Double I_ERI_G3xy_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[346];
  Double I_ERI_G3xz_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[347];
  Double I_ERI_G2x2y_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[348];
  Double I_ERI_G2xyz_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[349];
  Double I_ERI_G2x2z_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[350];
  Double I_ERI_Gx3y_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[351];
  Double I_ERI_Gx2yz_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[352];
  Double I_ERI_Gxy2z_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[353];
  Double I_ERI_Gx3z_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[354];
  Double I_ERI_G4y_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[355];
  Double I_ERI_G3yz_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[356];
  Double I_ERI_G2y2z_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[357];
  Double I_ERI_Gy3z_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[358];
  Double I_ERI_G4z_S_I4y2z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[359];
  Double I_ERI_G4x_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[360];
  Double I_ERI_G3xy_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[361];
  Double I_ERI_G3xz_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[362];
  Double I_ERI_G2x2y_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[363];
  Double I_ERI_G2xyz_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[364];
  Double I_ERI_G2x2z_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[365];
  Double I_ERI_Gx3y_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[366];
  Double I_ERI_Gx2yz_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[367];
  Double I_ERI_Gxy2z_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[368];
  Double I_ERI_Gx3z_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[369];
  Double I_ERI_G4y_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[370];
  Double I_ERI_G3yz_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[371];
  Double I_ERI_G2y2z_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[372];
  Double I_ERI_Gy3z_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[373];
  Double I_ERI_G4z_S_I3y3z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[374];
  Double I_ERI_G4x_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[375];
  Double I_ERI_G3xy_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[376];
  Double I_ERI_G3xz_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[377];
  Double I_ERI_G2x2y_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[378];
  Double I_ERI_G2xyz_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[379];
  Double I_ERI_G2x2z_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[380];
  Double I_ERI_Gx3y_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[381];
  Double I_ERI_Gx2yz_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[382];
  Double I_ERI_Gxy2z_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[383];
  Double I_ERI_Gx3z_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[384];
  Double I_ERI_G4y_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[385];
  Double I_ERI_G3yz_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[386];
  Double I_ERI_G2y2z_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[387];
  Double I_ERI_Gy3z_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[388];
  Double I_ERI_G4z_S_I2y4z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[389];
  Double I_ERI_G4x_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[390];
  Double I_ERI_G3xy_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[391];
  Double I_ERI_G3xz_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[392];
  Double I_ERI_G2x2y_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[393];
  Double I_ERI_G2xyz_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[394];
  Double I_ERI_G2x2z_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[395];
  Double I_ERI_Gx3y_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[396];
  Double I_ERI_Gx2yz_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[397];
  Double I_ERI_Gxy2z_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[398];
  Double I_ERI_Gx3z_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[399];
  Double I_ERI_G4y_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[400];
  Double I_ERI_G3yz_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[401];
  Double I_ERI_G2y2z_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[402];
  Double I_ERI_Gy3z_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[403];
  Double I_ERI_G4z_S_Iy5z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[404];
  Double I_ERI_G4x_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[405];
  Double I_ERI_G3xy_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[406];
  Double I_ERI_G3xz_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[407];
  Double I_ERI_G2x2y_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[408];
  Double I_ERI_G2xyz_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[409];
  Double I_ERI_G2x2z_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[410];
  Double I_ERI_Gx3y_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[411];
  Double I_ERI_Gx2yz_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[412];
  Double I_ERI_Gxy2z_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[413];
  Double I_ERI_Gx3z_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[414];
  Double I_ERI_G4y_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[415];
  Double I_ERI_G3yz_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[416];
  Double I_ERI_G2y2z_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[417];
  Double I_ERI_Gy3z_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[418];
  Double I_ERI_G4z_S_I6z_S_vrr = SQ_ERI_G_S_I_S_vrr_array[419];

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_I_S_C1003001003_cc
   * doing contraction work for VRR part 
   * totally 0 integrals are omitted 
   ************************************************************/
  Double SQ_ERI_G_S_I_S_C1003001003_cc_coefs = ic2_1*jc2_1*gamma*gamma;
  SQ_ERI_G_S_I_S_C1003001003_cc[0] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[1] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[2] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[3] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[4] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[5] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[6] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[7] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[8] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[9] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[10] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[11] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[12] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[13] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[14] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I6x_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[15] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[16] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[17] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[18] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[19] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[20] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[21] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[22] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[23] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[24] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[25] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[26] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[27] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[28] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[29] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I5xy_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[30] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[31] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[32] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[33] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[34] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[35] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[36] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[37] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[38] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[39] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[40] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[41] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[42] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[43] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[44] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I5xz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[45] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[46] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[47] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[48] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[49] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[50] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[51] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[52] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[53] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[54] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[55] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[56] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[57] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[58] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[59] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I4x2y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[60] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[61] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[62] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[63] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[64] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[65] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[66] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[67] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[68] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[69] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[70] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[71] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[72] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[73] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[74] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I4xyz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[75] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[76] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[77] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[78] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[79] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[80] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[81] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[82] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[83] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[84] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[85] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[86] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[87] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[88] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[89] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I4x2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[90] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[91] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[92] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[93] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[94] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[95] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[96] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[97] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[98] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[99] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[100] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[101] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[102] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[103] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[104] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I3x3y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[105] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[106] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[107] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[108] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[109] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[110] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[111] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[112] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[113] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[114] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[115] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[116] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[117] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[118] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[119] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I3x2yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[120] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[121] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[122] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[123] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[124] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[125] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[126] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[127] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[128] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[129] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[130] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[131] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[132] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[133] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[134] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I3xy2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[135] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[136] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[137] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[138] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[139] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[140] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[141] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[142] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[143] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[144] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[145] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[146] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[147] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[148] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[149] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I3x3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[150] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[151] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[152] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[153] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[154] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[155] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[156] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[157] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[158] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[159] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[160] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[161] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[162] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[163] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[164] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I2x4y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[165] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[166] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[167] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[168] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[169] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[170] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[171] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[172] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[173] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[174] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[175] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[176] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[177] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[178] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[179] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I2x3yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[180] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[181] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[182] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[183] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[184] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[185] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[186] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[187] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[188] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[189] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[190] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[191] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[192] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[193] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[194] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I2x2y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[195] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[196] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[197] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[198] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[199] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[200] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[201] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[202] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[203] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[204] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[205] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[206] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[207] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[208] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[209] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I2xy3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[210] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[211] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[212] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[213] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[214] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[215] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[216] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[217] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[218] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[219] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[220] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[221] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[222] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[223] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[224] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I2x4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[225] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[226] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[227] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[228] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[229] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[230] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[231] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[232] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[233] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[234] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[235] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[236] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[237] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[238] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[239] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_Ix5y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[240] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[241] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[242] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[243] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[244] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[245] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[246] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[247] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[248] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[249] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[250] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[251] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[252] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[253] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[254] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_Ix4yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[255] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[256] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[257] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[258] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[259] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[260] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[261] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[262] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[263] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[264] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[265] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[266] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[267] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[268] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[269] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_Ix3y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[270] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[271] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[272] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[273] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[274] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[275] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[276] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[277] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[278] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[279] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[280] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[281] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[282] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[283] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[284] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_Ix2y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[285] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[286] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[287] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[288] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[289] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[290] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[291] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[292] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[293] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[294] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[295] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[296] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[297] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[298] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[299] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_Ixy4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[300] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[301] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[302] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[303] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[304] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[305] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[306] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[307] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[308] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[309] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[310] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[311] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[312] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[313] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[314] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_Ix5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[315] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[316] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[317] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[318] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[319] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[320] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[321] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[322] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[323] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[324] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[325] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[326] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[327] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[328] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[329] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I6y_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[330] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[331] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[332] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[333] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[334] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[335] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[336] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[337] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[338] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[339] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[340] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[341] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[342] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[343] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[344] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I5yz_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[345] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[346] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[347] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[348] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[349] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[350] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[351] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[352] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[353] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[354] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[355] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[356] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[357] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[358] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[359] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I4y2z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[360] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[361] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[362] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[363] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[364] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[365] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[366] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[367] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[368] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[369] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[370] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[371] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[372] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[373] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[374] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I3y3z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[375] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[376] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[377] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[378] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[379] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[380] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[381] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[382] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[383] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[384] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[385] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[386] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[387] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[388] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[389] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I2y4z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[390] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[391] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[392] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[393] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[394] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[395] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[396] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[397] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[398] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[399] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[400] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[401] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[402] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[403] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[404] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_Iy5z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[405] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4x_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[406] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xy_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[407] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3xz_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[408] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2y_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[409] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2xyz_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[410] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2x2z_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[411] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3y_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[412] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx2yz_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[413] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gxy2z_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[414] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gx3z_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[415] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4y_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[416] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G3yz_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[417] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G2y2z_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[418] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_Gy3z_S_I6z_S_vrr;
  SQ_ERI_G_S_I_S_C1003001003_cc[419] += SQ_ERI_G_S_I_S_C1003001003_cc_coefs*I_ERI_G4z_S_I6z_S_vrr;
}
