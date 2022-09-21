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
#include <boost/math/special_functions/gamma.hpp>
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
// BRA1 as redundant position, total RHS integrals evaluated as: 1172654
// BRA2 as redundant position, total RHS integrals evaluated as: 1127042
// KET1 as redundant position, total RHS integrals evaluated as: 1172461
// KET2 as redundant position, total RHS integrals evaluated as: 1127042
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


// VRR working function: 1
void hgp_os_eri_f_d_f_d_d2_vrr_1(const Double& oned2k, const Double& oned2z, const Double& rhod2zsq, const Double& WPX, const Double& WPY, const Double& WPZ, const Double& oned2e, const Double& rhod2esq, const Double& WQX, const Double& WQY, const Double& WQZ, const Double& PAX, const Double& PAY, const Double& PAZ, const Double& QCX, const Double& QCY, const Double& QCZ, const Double& alpha, const Double& beta, const Double& gamma, const Double& delta, const Double& I_ERI_S_S_S_S_vrr, const Double& I_ERI_S_S_S_S_M1_vrr, const Double& I_ERI_S_S_S_S_M2_vrr, const Double& I_ERI_S_S_S_S_M3_vrr, const Double& I_ERI_S_S_S_S_M4_vrr, const Double& I_ERI_S_S_S_S_M5_vrr, const Double& I_ERI_S_S_S_S_M6_vrr, const Double& I_ERI_S_S_S_S_M7_vrr, const Double& I_ERI_S_S_S_S_M8_vrr, const Double& I_ERI_S_S_S_S_M9_vrr, const Double& I_ERI_S_S_S_S_M10_vrr, const Double& I_ERI_S_S_S_S_M11_vrr, const Double& I_ERI_S_S_S_S_M12_vrr, Double* SQ_ERI_D_S_D_S_vrr_array, Double* SQ_ERI_P_S_F_S_vrr_array, Double* SQ_ERI_F_S_P_S_vrr_array, Double* SQ_ERI_D_S_F_S_vrr_array, Double* SQ_ERI_F_S_D_S_vrr_array, Double* SQ_ERI_P_S_G_S_vrr_array, Double* SQ_ERI_G_S_P_S_vrr_array, Double* SQ_ERI_F_S_F_S_vrr_array, Double* SQ_ERI_D_S_G_S_vrr_array, Double* SQ_ERI_G_S_D_S_vrr_array, Double* SQ_ERI_P_S_H_S_vrr_array, Double* SQ_ERI_H_S_P_S_vrr_array, Double* SQ_ERI_F_S_G_S_vrr_array, Double* SQ_ERI_G_S_F_S_vrr_array, Double* SQ_ERI_D_S_H_S_vrr_array, Double* SQ_ERI_H_S_D_S_vrr_array, Double* SQ_ERI_G_S_G_S_vrr_array, Double* SQ_ERI_F_S_H_S_vrr_array, Double* SQ_ERI_H_S_F_S_vrr_array, Double* SQ_ERI_D_S_I_S_vrr_array, Double* SQ_ERI_I_S_D_S_vrr_array, Double* SQ_ERI_G_S_H_S_vrr_array, Double* SQ_ERI_H_S_G_S_vrr_array, Double* SQ_ERI_F_S_I_S_vrr_array, Double* SQ_ERI_I_S_F_S_vrr_array, Double* SQ_ERI_H_S_H_S_vrr_array, Double* SQ_ERI_G_S_I_S_vrr_array, Double* SQ_ERI_I_S_G_S_vrr_array, Double* SQ_ERI_F_S_K_S_vrr_array, Double* SQ_ERI_K_S_F_S_vrr_array, Double* SQ_ERI_H_S_I_S_vrr_array, Double* SQ_ERI_I_S_H_S_vrr_array, Double* SQ_ERI_G_S_K_S_vrr_array, Double* SQ_ERI_K_S_G_S_vrr_array, Double* SQ_ERI_I_S_I_S_vrr_array, Double* SQ_ERI_H_S_K_S_vrr_array, Double* SQ_ERI_K_S_H_S_vrr_array  );


// HRR working function: 1
void hgp_os_eri_f_d_f_d_d2_hrr1_1(const Double* C, const Double* D, const Double* SQ_ERI_F_S_D_S, const Double* SQ_ERI_F_S_P_S, const Double* SQ_ERI_G_S_D_S, const Double* SQ_ERI_G_S_P_S, const Double* SQ_ERI_H_S_D_S, const Double* SQ_ERI_H_S_P_S, const Double* SQ_ERI_D_S_F_S, const Double* SQ_ERI_D_S_D_S, const Double* SQ_ERI_F_S_F_S, const Double* SQ_ERI_G_S_F_S, const Double* SQ_ERI_H_S_F_S, const Double* SQ_ERI_P_S_G_S, const Double* SQ_ERI_P_S_F_S, const Double* SQ_ERI_D_S_G_S, const Double* SQ_ERI_F_S_G_S, const Double* SQ_ERI_G_S_G_S, const Double* SQ_ERI_P_S_H_S, const Double* SQ_ERI_D_S_H_S, const Double* SQ_ERI_F_S_H_S, const Double* SQ_ERI_G_S_F_S_a, const Double* SQ_ERI_G_S_D_S_a, const Double* SQ_ERI_H_S_F_S_a, const Double* SQ_ERI_H_S_D_S_a, const Double* SQ_ERI_I_S_F_S_a, const Double* SQ_ERI_I_S_D_S_a, const Double* SQ_ERI_F_S_G_S_a, const Double* SQ_ERI_F_S_F_S_a, const Double* SQ_ERI_G_S_G_S_a, const Double* SQ_ERI_H_S_G_S_a, const Double* SQ_ERI_I_S_G_S_a, const Double* SQ_ERI_F_S_H_S_a, const Double* SQ_ERI_G_S_H_S_a, const Double* SQ_ERI_H_S_H_S_a, const Double* SQ_ERI_F_S_F_S_b, const Double* SQ_ERI_F_S_D_S_b, const Double* SQ_ERI_G_S_F_S_b, const Double* SQ_ERI_G_S_D_S_b, const Double* SQ_ERI_H_S_F_S_b, const Double* SQ_ERI_H_S_D_S_b, const Double* SQ_ERI_I_S_F_S_b, const Double* SQ_ERI_I_S_D_S_b, const Double* SQ_ERI_D_S_G_S_b, const Double* SQ_ERI_D_S_F_S_b, const Double* SQ_ERI_F_S_G_S_b, const Double* SQ_ERI_G_S_G_S_b, const Double* SQ_ERI_H_S_G_S_b, const Double* SQ_ERI_I_S_G_S_b, const Double* SQ_ERI_D_S_H_S_b, const Double* SQ_ERI_F_S_H_S_b, const Double* SQ_ERI_G_S_H_S_b, const Double* SQ_ERI_H_S_H_S_b, const Double* SQ_ERI_F_S_G_S_c, const Double* SQ_ERI_F_S_F_S_c, const Double* SQ_ERI_G_S_G_S_c, const Double* SQ_ERI_G_S_F_S_c, const Double* SQ_ERI_H_S_G_S_c, const Double* SQ_ERI_H_S_F_S_c, const Double* SQ_ERI_D_S_H_S_c, const Double* SQ_ERI_D_S_G_S_c, const Double* SQ_ERI_F_S_H_S_c, const Double* SQ_ERI_G_S_H_S_c, const Double* SQ_ERI_H_S_H_S_c, const Double* SQ_ERI_D_S_I_S_c, const Double* SQ_ERI_F_S_I_S_c, const Double* SQ_ERI_G_S_I_S_c, const Double* SQ_ERI_H_S_G_S_aa, const Double* SQ_ERI_H_S_F_S_aa, const Double* SQ_ERI_I_S_G_S_aa, const Double* SQ_ERI_I_S_F_S_aa, const Double* SQ_ERI_K_S_G_S_aa, const Double* SQ_ERI_K_S_F_S_aa, const Double* SQ_ERI_H_S_H_S_aa, const Double* SQ_ERI_I_S_H_S_aa, const Double* SQ_ERI_K_S_H_S_aa, const Double* SQ_ERI_G_S_G_S_ab, const Double* SQ_ERI_G_S_F_S_ab, const Double* SQ_ERI_H_S_G_S_ab, const Double* SQ_ERI_H_S_F_S_ab, Double* SQ_ERI_F_S_P_D, Double* SQ_ERI_G_S_P_D, Double* SQ_ERI_H_S_P_D, Double* SQ_ERI_D_S_D_D, Double* SQ_ERI_F_S_D_D, Double* SQ_ERI_G_S_D_D, Double* SQ_ERI_P_S_F_D, Double* SQ_ERI_D_S_F_D, Double* SQ_ERI_F_S_F_D, Double* SQ_ERI_G_S_D_D_a, Double* SQ_ERI_H_S_D_D_a, Double* SQ_ERI_I_S_D_D_a, Double* SQ_ERI_F_S_F_D_a, Double* SQ_ERI_G_S_F_D_a, Double* SQ_ERI_H_S_F_D_a, Double* SQ_ERI_F_S_D_D_b, Double* SQ_ERI_G_S_D_D_b, Double* SQ_ERI_H_S_D_D_b, Double* SQ_ERI_I_S_D_D_b, Double* SQ_ERI_D_S_F_D_b, Double* SQ_ERI_F_S_F_D_b, Double* SQ_ERI_G_S_F_D_b, Double* SQ_ERI_H_S_F_D_b, Double* SQ_ERI_F_S_F_D_c, Double* SQ_ERI_G_S_F_D_c, Double* SQ_ERI_H_S_F_D_c, Double* SQ_ERI_D_S_G_D_c, Double* SQ_ERI_F_S_G_D_c, Double* SQ_ERI_G_S_G_D_c, Double* SQ_ERI_H_S_F_D_aa, Double* SQ_ERI_I_S_F_D_aa, Double* SQ_ERI_K_S_F_D_aa, Double* SQ_ERI_G_S_F_P_ab, Double* SQ_ERI_H_S_F_P_ab  );


// HRR working function: 2
void hgp_os_eri_f_d_f_d_d2_hrr1_2(const Double* C, const Double* D, const Double* SQ_ERI_I_S_G_S_ab, const Double* SQ_ERI_I_S_F_S_ab, const Double* SQ_ERI_K_S_G_S_ab, const Double* SQ_ERI_K_S_F_S_ab, const Double* SQ_ERI_G_S_H_S_ab, const Double* SQ_ERI_G_S_G_S_ab, const Double* SQ_ERI_H_S_H_S_ab, const Double* SQ_ERI_H_S_G_S_ab, const Double* SQ_ERI_I_S_H_S_ab, const Double* SQ_ERI_K_S_H_S_ab, const Double* SQ_ERI_G_S_F_P_ab, const Double* SQ_ERI_H_S_F_P_ab, const Double* SQ_ERI_G_S_H_S_ac, const Double* SQ_ERI_G_S_G_S_ac, const Double* SQ_ERI_H_S_H_S_ac, const Double* SQ_ERI_H_S_G_S_ac, const Double* SQ_ERI_I_S_H_S_ac, const Double* SQ_ERI_I_S_G_S_ac, const Double* SQ_ERI_G_S_I_S_ac, const Double* SQ_ERI_H_S_I_S_ac, const Double* SQ_ERI_I_S_I_S_ac, const Double* SQ_ERI_F_S_G_S_bb, const Double* SQ_ERI_F_S_F_S_bb, const Double* SQ_ERI_G_S_G_S_bb, const Double* SQ_ERI_G_S_F_S_bb, const Double* SQ_ERI_H_S_G_S_bb, const Double* SQ_ERI_H_S_F_S_bb, const Double* SQ_ERI_I_S_G_S_bb, const Double* SQ_ERI_I_S_F_S_bb, const Double* SQ_ERI_K_S_G_S_bb, const Double* SQ_ERI_K_S_F_S_bb, const Double* SQ_ERI_F_S_H_S_bb, const Double* SQ_ERI_G_S_H_S_bb, const Double* SQ_ERI_H_S_H_S_bb, const Double* SQ_ERI_I_S_H_S_bb, const Double* SQ_ERI_K_S_H_S_bb, const Double* SQ_ERI_F_S_H_S_bc, const Double* SQ_ERI_F_S_G_S_bc, const Double* SQ_ERI_G_S_H_S_bc, const Double* SQ_ERI_G_S_G_S_bc, const Double* SQ_ERI_H_S_H_S_bc, const Double* SQ_ERI_H_S_G_S_bc, const Double* SQ_ERI_I_S_H_S_bc, const Double* SQ_ERI_I_S_G_S_bc, const Double* SQ_ERI_F_S_I_S_bc, const Double* SQ_ERI_G_S_I_S_bc, const Double* SQ_ERI_H_S_I_S_bc, const Double* SQ_ERI_I_S_I_S_bc, Double* SQ_ERI_G_S_F_D_ab, Double* SQ_ERI_H_S_F_D_ab, Double* SQ_ERI_I_S_F_D_ab, Double* SQ_ERI_K_S_F_D_ab, Double* SQ_ERI_G_S_G_D_ac, Double* SQ_ERI_H_S_G_D_ac, Double* SQ_ERI_I_S_G_D_ac, Double* SQ_ERI_F_S_F_D_bb, Double* SQ_ERI_G_S_F_D_bb, Double* SQ_ERI_H_S_F_D_bb, Double* SQ_ERI_I_S_F_D_bb, Double* SQ_ERI_K_S_F_D_bb, Double* SQ_ERI_F_S_G_D_bc, Double* SQ_ERI_G_S_G_D_bc, Double* SQ_ERI_H_S_G_D_bc, Double* SQ_ERI_I_S_G_D_bc  );


// HRR working function: 3
void hgp_os_eri_f_d_f_d_d2_hrr1_3(const Double* C, const Double* D, const Double* SQ_ERI_F_S_I_S_cc, const Double* SQ_ERI_F_S_H_S_cc, const Double* SQ_ERI_G_S_I_S_cc, const Double* SQ_ERI_G_S_H_S_cc, const Double* SQ_ERI_H_S_I_S_cc, const Double* SQ_ERI_H_S_H_S_cc, const Double* SQ_ERI_F_S_K_S_cc, const Double* SQ_ERI_G_S_K_S_cc, const Double* SQ_ERI_H_S_K_S_cc, Double* SQ_ERI_F_S_H_D_cc, Double* SQ_ERI_G_S_H_D_cc, Double* SQ_ERI_H_S_H_D_cc  );


// HRR working function: 1
void hgp_os_eri_f_d_f_d_d2_hrr2_1(const Double* A, const Double* B, const Double* SQ_ERI_D_S_F_D, const Double* SQ_ERI_P_S_F_D, const Double* SQ_ERI_F_S_D_D, const Double* SQ_ERI_D_S_D_D, const Double* SQ_ERI_F_S_F_D, const Double* SQ_ERI_G_S_P_D, const Double* SQ_ERI_F_S_P_D, const Double* SQ_ERI_G_S_D_D, const Double* SQ_ERI_H_S_P_D, const Double* SQ_ERI_G_S_F_D_a, const Double* SQ_ERI_F_S_F_D_a, const Double* SQ_ERI_H_S_D_D_a, const Double* SQ_ERI_G_S_D_D_a, const Double* SQ_ERI_H_S_F_D_a, const Double* SQ_ERI_I_S_D_D_a, const Double* SQ_ERI_F_S_F_D_b, const Double* SQ_ERI_D_S_F_D_b, const Double* SQ_ERI_G_S_D_D_b, const Double* SQ_ERI_F_S_D_D_b, const Double* SQ_ERI_G_S_F_D_b, const Double* SQ_ERI_H_S_D_D_b, const Double* SQ_ERI_H_S_F_D_b, const Double* SQ_ERI_I_S_D_D_b, const Double* SQ_ERI_F_S_G_D_c, const Double* SQ_ERI_D_S_G_D_c, const Double* SQ_ERI_G_S_F_D_c, const Double* SQ_ERI_F_S_F_D_c, Double* SQ_ERI_D_P_F_D, Double* SQ_ERI_P_D_F_D, Double* SQ_ERI_F_P_D_D, Double* SQ_ERI_D_D_D_D, Double* SQ_ERI_F_D_P_D, Double* SQ_ERI_G_P_F_D_a, Double* SQ_ERI_F_D_F_D_a, Double* SQ_ERI_G_D_D_D_a, Double* SQ_ERI_F_D_F_D_b, Double* SQ_ERI_D_F_F_D_b, Double* SQ_ERI_F_F_D_D_b, Double* SQ_ERI_D_P_G_D_c, Double* SQ_ERI_F_P_F_D_c  );


// HRR working function: 2
void hgp_os_eri_f_d_f_d_d2_hrr2_2(const Double* A, const Double* B, const Double* SQ_ERI_G_S_G_D_c, const Double* SQ_ERI_F_S_G_D_c, const Double* SQ_ERI_D_P_G_D_c, const Double* SQ_ERI_H_S_F_D_c, const Double* SQ_ERI_G_S_F_D_c, const Double* SQ_ERI_F_P_F_D_c, const Double* SQ_ERI_I_S_F_D_aa, const Double* SQ_ERI_H_S_F_D_aa, const Double* SQ_ERI_K_S_F_D_aa, const Double* SQ_ERI_H_S_F_D_ab, const Double* SQ_ERI_G_S_F_D_ab, const Double* SQ_ERI_I_S_F_D_ab, const Double* SQ_ERI_K_S_F_D_ab, Double* SQ_ERI_F_P_G_D_c, Double* SQ_ERI_D_D_G_D_c, Double* SQ_ERI_F_D_F_D_c, Double* SQ_ERI_H_D_F_D_aa, Double* SQ_ERI_G_F_F_D_ab  );


// HRR working function: 3
void hgp_os_eri_f_d_f_d_d2_hrr2_3(const Double* A, const Double* B, const Double* SQ_ERI_H_S_G_D_ac, const Double* SQ_ERI_G_S_G_D_ac, const Double* SQ_ERI_I_S_G_D_ac, const Double* SQ_ERI_G_S_F_D_bb, const Double* SQ_ERI_F_S_F_D_bb, const Double* SQ_ERI_H_S_F_D_bb, const Double* SQ_ERI_I_S_F_D_bb, const Double* SQ_ERI_K_S_F_D_bb, Double* SQ_ERI_G_D_G_D_ac, Double* SQ_ERI_F_G_F_D_bb  );


// HRR working function: 4
void hgp_os_eri_f_d_f_d_d2_hrr2_4(const Double* A, const Double* B, const Double* SQ_ERI_G_S_G_D_bc, const Double* SQ_ERI_F_S_G_D_bc, const Double* SQ_ERI_H_S_G_D_bc, const Double* SQ_ERI_I_S_G_D_bc, const Double* SQ_ERI_G_S_H_D_cc, const Double* SQ_ERI_F_S_H_D_cc, const Double* SQ_ERI_H_S_H_D_cc, Double* SQ_ERI_F_F_G_D_bc, Double* SQ_ERI_F_D_H_D_cc  );


// NON-RR working function: 1
void hgp_os_eri_f_d_f_d_d2_deriv_1(const Double* SQ_ERI_H_D_F_D_aa, const Double* SQ_ERI_F_D_F_D_a, const Double* SQ_ERI_P_D_F_D, const Double* SQ_ERI_G_F_F_D_ab, const Double* SQ_ERI_G_P_F_D_a, const Double* SQ_ERI_D_F_F_D_b, const Double* SQ_ERI_D_P_F_D, Double* abcd);


// NON-RR working function: 2
void hgp_os_eri_f_d_f_d_d2_deriv_2(const Double* SQ_ERI_G_F_F_D_ab, const Double* SQ_ERI_G_P_F_D_a, const Double* SQ_ERI_D_F_F_D_b, const Double* SQ_ERI_D_P_F_D, const Double* SQ_ERI_G_D_G_D_ac, const Double* SQ_ERI_G_D_D_D_a, const Double* SQ_ERI_D_D_G_D_c, const Double* SQ_ERI_D_D_D_D, Double* abcd);


// NON-RR working function: 3
void hgp_os_eri_f_d_f_d_d2_deriv_3(const Double* SQ_ERI_F_G_F_D_bb, const Double* SQ_ERI_F_D_F_D_b, const Double* SQ_ERI_F_S_F_D, const Double* SQ_ERI_F_F_G_D_bc, const Double* SQ_ERI_F_F_D_D_b, const Double* SQ_ERI_F_P_G_D_c, const Double* SQ_ERI_F_P_D_D, Double* abcd);


// NON-RR working function: 4
void hgp_os_eri_f_d_f_d_d2_deriv_4(const Double* SQ_ERI_F_F_G_D_bc, const Double* SQ_ERI_F_F_D_D_b, const Double* SQ_ERI_F_P_G_D_c, const Double* SQ_ERI_F_P_D_D, const Double* SQ_ERI_F_D_H_D_cc, const Double* SQ_ERI_F_D_F_D_c, const Double* SQ_ERI_F_D_P_D, Double* abcd);


// VRR contraction function: 1
void hgp_os_eri_f_d_f_d_d2_vrrcont_1(const Double& alpha, const Double& beta, const Double& gamma, const Double& delta, const Double* SQ_ERI_F_S_H_S_vrr_array, const Double* SQ_ERI_F_S_G_S_vrr_array, const Double* SQ_ERI_F_S_F_S_vrr_array, const Double* SQ_ERI_H_S_H_S_vrr_array, const Double* SQ_ERI_H_S_G_S_vrr_array, const Double* SQ_ERI_H_S_F_S_vrr_array, const Double* SQ_ERI_G_S_H_S_vrr_array, const Double* SQ_ERI_G_S_G_S_vrr_array, Double* SQ_ERI_F_S_H_S, Double* SQ_ERI_F_S_G_S, Double* SQ_ERI_F_S_F_S, Double* SQ_ERI_H_S_H_S_a, Double* SQ_ERI_H_S_G_S_a, Double* SQ_ERI_H_S_F_S_a, Double* SQ_ERI_G_S_H_S_a, Double* SQ_ERI_G_S_G_S_a, Double* SQ_ERI_G_S_H_S_c, Double* SQ_ERI_G_S_G_S_c, Double* SQ_ERI_F_S_H_S_c, Double* SQ_ERI_F_S_G_S_c, Double* SQ_ERI_G_S_G_S, Double* SQ_ERI_H_S_H_S_aa, Double* SQ_ERI_H_S_G_S_aa, Double* SQ_ERI_H_S_F_S_aa, Double* SQ_ERI_F_S_H_S_a, Double* SQ_ERI_F_S_G_S_a, Double* SQ_ERI_F_S_F_S_a, Double* SQ_ERI_H_S_H_S_ac, Double* SQ_ERI_H_S_G_S_ac, Double* SQ_ERI_G_S_H_S_ac, Double* SQ_ERI_G_S_G_S_ac, Double* SQ_ERI_H_S_H_S_b, Double* SQ_ERI_H_S_G_S_b, Double* SQ_ERI_H_S_F_S_b, Double* SQ_ERI_G_S_H_S_b, Double* SQ_ERI_G_S_G_S_b, Double* SQ_ERI_F_S_H_S_b, Double* SQ_ERI_F_S_G_S_b, Double* SQ_ERI_F_S_F_S_b, Double* SQ_ERI_H_S_H_S_cc, Double* SQ_ERI_G_S_H_S_cc, Double* SQ_ERI_F_S_H_S_cc, Double* SQ_ERI_H_S_H_S_c, Double* SQ_ERI_H_S_G_S_c, Double* SQ_ERI_H_S_F_S_c, Double* SQ_ERI_F_S_F_S_c, Double* SQ_ERI_H_S_F_S, Double* SQ_ERI_H_S_H_S_ab, Double* SQ_ERI_H_S_G_S_ab, Double* SQ_ERI_H_S_F_S_ab, Double* SQ_ERI_G_S_H_S_ab, Double* SQ_ERI_G_S_G_S_ab, Double* SQ_ERI_H_S_H_S_bc, Double* SQ_ERI_H_S_G_S_bc, Double* SQ_ERI_G_S_H_S_bc, Double* SQ_ERI_G_S_G_S_bc, Double* SQ_ERI_F_S_H_S_bc, Double* SQ_ERI_F_S_G_S_bc, Double* SQ_ERI_H_S_H_S_bb, Double* SQ_ERI_H_S_G_S_bb, Double* SQ_ERI_H_S_F_S_bb, Double* SQ_ERI_G_S_H_S_bb, Double* SQ_ERI_G_S_G_S_bb, Double* SQ_ERI_F_S_H_S_bb, Double* SQ_ERI_F_S_G_S_bb, Double* SQ_ERI_F_S_F_S_bb  );


// VRR contraction function: 2
void hgp_os_eri_f_d_f_d_d2_vrrcont_2(const Double& alpha, const Double& beta, const Double& gamma, const Double& delta, const Double* SQ_ERI_G_S_F_S_vrr_array, const Double* SQ_ERI_D_S_H_S_vrr_array, const Double* SQ_ERI_D_S_G_S_vrr_array, const Double* SQ_ERI_D_S_F_S_vrr_array, const Double* SQ_ERI_G_S_I_S_vrr_array, const Double* SQ_ERI_F_S_I_S_vrr_array, const Double* SQ_ERI_G_S_D_S_vrr_array, const Double* SQ_ERI_F_S_D_S_vrr_array, const Double* SQ_ERI_K_S_H_S_vrr_array, const Double* SQ_ERI_K_S_G_S_vrr_array, const Double* SQ_ERI_K_S_F_S_vrr_array, const Double* SQ_ERI_I_S_H_S_vrr_array, const Double* SQ_ERI_I_S_G_S_vrr_array, const Double* SQ_ERI_I_S_F_S_vrr_array, Double* SQ_ERI_G_S_F_S_a, Double* SQ_ERI_D_S_H_S, Double* SQ_ERI_D_S_G_S, Double* SQ_ERI_D_S_F_S, Double* SQ_ERI_G_S_I_S_c, Double* SQ_ERI_F_S_I_S_c, Double* SQ_ERI_G_S_F_S, Double* SQ_ERI_G_S_D_S, Double* SQ_ERI_F_S_D_S, Double* SQ_ERI_K_S_H_S_aa, Double* SQ_ERI_K_S_G_S_aa, Double* SQ_ERI_K_S_F_S_aa, Double* SQ_ERI_I_S_H_S_aa, Double* SQ_ERI_I_S_G_S_aa, Double* SQ_ERI_I_S_F_S_aa, Double* SQ_ERI_I_S_H_S_ac, Double* SQ_ERI_I_S_G_S_ac, Double* SQ_ERI_G_S_I_S_ac, Double* SQ_ERI_I_S_G_S_a, Double* SQ_ERI_I_S_F_S_a, Double* SQ_ERI_G_S_D_S_a, Double* SQ_ERI_D_S_H_S_c, Double* SQ_ERI_D_S_G_S_c, Double* SQ_ERI_G_S_F_S_b, Double* SQ_ERI_G_S_I_S_cc, Double* SQ_ERI_F_S_I_S_cc, Double* SQ_ERI_G_S_F_S_c, Double* SQ_ERI_K_S_H_S_ab, Double* SQ_ERI_K_S_G_S_ab, Double* SQ_ERI_K_S_F_S_ab, Double* SQ_ERI_I_S_H_S_ab, Double* SQ_ERI_I_S_G_S_ab, Double* SQ_ERI_I_S_F_S_ab, Double* SQ_ERI_G_S_F_S_ab, Double* SQ_ERI_D_S_H_S_b, Double* SQ_ERI_D_S_G_S_b, Double* SQ_ERI_D_S_F_S_b, Double* SQ_ERI_I_S_H_S_bc, Double* SQ_ERI_I_S_G_S_bc, Double* SQ_ERI_G_S_I_S_bc, Double* SQ_ERI_F_S_I_S_bc, Double* SQ_ERI_I_S_G_S_b, Double* SQ_ERI_I_S_F_S_b, Double* SQ_ERI_G_S_D_S_b, Double* SQ_ERI_F_S_D_S_b, Double* SQ_ERI_K_S_H_S_bb, Double* SQ_ERI_K_S_G_S_bb, Double* SQ_ERI_K_S_F_S_bb, Double* SQ_ERI_I_S_H_S_bb, Double* SQ_ERI_I_S_G_S_bb, Double* SQ_ERI_I_S_F_S_bb, Double* SQ_ERI_G_S_F_S_bb  );


// VRR contraction function: 3
void hgp_os_eri_f_d_f_d_d2_vrrcont_3(const Double& alpha, const Double& beta, const Double& gamma, const Double& delta, const Double* SQ_ERI_P_S_H_S_vrr_array, const Double* SQ_ERI_P_S_G_S_vrr_array, const Double* SQ_ERI_P_S_F_S_vrr_array, const Double* SQ_ERI_I_S_I_S_vrr_array, const Double* SQ_ERI_H_S_I_S_vrr_array, const Double* SQ_ERI_I_S_D_S_vrr_array, const Double* SQ_ERI_H_S_D_S_vrr_array, const Double* SQ_ERI_D_S_I_S_vrr_array, const Double* SQ_ERI_D_S_D_S_vrr_array, const Double* SQ_ERI_H_S_K_S_vrr_array, const Double* SQ_ERI_G_S_K_S_vrr_array, const Double* SQ_ERI_F_S_K_S_vrr_array, const Double* SQ_ERI_H_S_P_S_vrr_array, const Double* SQ_ERI_G_S_P_S_vrr_array, const Double* SQ_ERI_F_S_P_S_vrr_array, Double* SQ_ERI_P_S_H_S, Double* SQ_ERI_P_S_G_S, Double* SQ_ERI_P_S_F_S, Double* SQ_ERI_I_S_I_S_ac, Double* SQ_ERI_H_S_I_S_ac, Double* SQ_ERI_I_S_D_S_a, Double* SQ_ERI_H_S_D_S_a, Double* SQ_ERI_D_S_I_S_c, Double* SQ_ERI_D_S_D_S, Double* SQ_ERI_H_S_K_S_cc, Double* SQ_ERI_H_S_I_S_cc, Double* SQ_ERI_G_S_K_S_cc, Double* SQ_ERI_F_S_K_S_cc, Double* SQ_ERI_H_S_D_S, Double* SQ_ERI_H_S_P_S, Double* SQ_ERI_G_S_P_S, Double* SQ_ERI_F_S_P_S, Double* SQ_ERI_I_S_I_S_bc, Double* SQ_ERI_H_S_I_S_bc, Double* SQ_ERI_I_S_D_S_b, Double* SQ_ERI_H_S_D_S_b  );

void hgp_os_eri_f_d_f_d_d2(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd, LocalMemScr& scr)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double* SQ_ERI_F_S_H_S = scr.getNewMemPos(210);
  Double* SQ_ERI_F_S_G_S = scr.getNewMemPos(150);
  Double* SQ_ERI_F_S_F_S = scr.getNewMemPos(100);
  Double* SQ_ERI_H_S_H_S_a = scr.getNewMemPos(441);
  Double* SQ_ERI_H_S_G_S_a = scr.getNewMemPos(315);
  Double* SQ_ERI_H_S_F_S_a = scr.getNewMemPos(210);
  Double* SQ_ERI_G_S_H_S_a = scr.getNewMemPos(315);
  Double* SQ_ERI_G_S_G_S_a = scr.getNewMemPos(225);
  Double* SQ_ERI_G_S_F_S_a = scr.getNewMemPos(150);
  Double* SQ_ERI_D_S_H_S = scr.getNewMemPos(126);
  Double* SQ_ERI_D_S_G_S = scr.getNewMemPos(90);
  Double* SQ_ERI_D_S_F_S = scr.getNewMemPos(60);
  Double* SQ_ERI_G_S_I_S_c = scr.getNewMemPos(420);
  Double* SQ_ERI_G_S_H_S_c = scr.getNewMemPos(315);
  Double* SQ_ERI_G_S_G_S_c = scr.getNewMemPos(225);
  Double* SQ_ERI_F_S_I_S_c = scr.getNewMemPos(280);
  Double* SQ_ERI_F_S_H_S_c = scr.getNewMemPos(210);
  Double* SQ_ERI_F_S_G_S_c = scr.getNewMemPos(150);
  Double* SQ_ERI_G_S_G_S = scr.getNewMemPos(225);
  Double* SQ_ERI_G_S_F_S = scr.getNewMemPos(150);
  Double* SQ_ERI_G_S_D_S = scr.getNewMemPos(90);
  Double* SQ_ERI_F_S_D_S = scr.getNewMemPos(60);
  Double* SQ_ERI_K_S_H_S_aa = scr.getNewMemPos(756);
  Double* SQ_ERI_K_S_G_S_aa = scr.getNewMemPos(540);
  Double* SQ_ERI_K_S_F_S_aa = scr.getNewMemPos(360);
  Double* SQ_ERI_I_S_H_S_aa = scr.getNewMemPos(588);
  Double* SQ_ERI_I_S_G_S_aa = scr.getNewMemPos(420);
  Double* SQ_ERI_I_S_F_S_aa = scr.getNewMemPos(280);
  Double* SQ_ERI_H_S_H_S_aa = scr.getNewMemPos(441);
  Double* SQ_ERI_H_S_G_S_aa = scr.getNewMemPos(315);
  Double* SQ_ERI_H_S_F_S_aa = scr.getNewMemPos(210);
  Double* SQ_ERI_F_S_H_S_a = scr.getNewMemPos(210);
  Double* SQ_ERI_F_S_G_S_a = scr.getNewMemPos(150);
  Double* SQ_ERI_F_S_F_S_a = scr.getNewMemPos(100);
  Double* SQ_ERI_P_S_H_S = scr.getNewMemPos(63);
  Double* SQ_ERI_P_S_G_S = scr.getNewMemPos(45);
  Double* SQ_ERI_P_S_F_S = scr.getNewMemPos(30);
  Double* SQ_ERI_I_S_I_S_ac = scr.getNewMemPos(784);
  Double* SQ_ERI_I_S_H_S_ac = scr.getNewMemPos(588);
  Double* SQ_ERI_I_S_G_S_ac = scr.getNewMemPos(420);
  Double* SQ_ERI_H_S_I_S_ac = scr.getNewMemPos(588);
  Double* SQ_ERI_H_S_H_S_ac = scr.getNewMemPos(441);
  Double* SQ_ERI_H_S_G_S_ac = scr.getNewMemPos(315);
  Double* SQ_ERI_G_S_I_S_ac = scr.getNewMemPos(420);
  Double* SQ_ERI_G_S_H_S_ac = scr.getNewMemPos(315);
  Double* SQ_ERI_G_S_G_S_ac = scr.getNewMemPos(225);
  Double* SQ_ERI_I_S_G_S_a = scr.getNewMemPos(420);
  Double* SQ_ERI_I_S_F_S_a = scr.getNewMemPos(280);
  Double* SQ_ERI_I_S_D_S_a = scr.getNewMemPos(168);
  Double* SQ_ERI_H_S_D_S_a = scr.getNewMemPos(126);
  Double* SQ_ERI_G_S_D_S_a = scr.getNewMemPos(90);
  Double* SQ_ERI_D_S_I_S_c = scr.getNewMemPos(168);
  Double* SQ_ERI_D_S_H_S_c = scr.getNewMemPos(126);
  Double* SQ_ERI_D_S_G_S_c = scr.getNewMemPos(90);
  Double* SQ_ERI_D_S_D_S = scr.getNewMemPos(36);
  Double* SQ_ERI_H_S_H_S_b = scr.getNewMemPos(441);
  Double* SQ_ERI_H_S_G_S_b = scr.getNewMemPos(315);
  Double* SQ_ERI_H_S_F_S_b = scr.getNewMemPos(210);
  Double* SQ_ERI_G_S_H_S_b = scr.getNewMemPos(315);
  Double* SQ_ERI_G_S_G_S_b = scr.getNewMemPos(225);
  Double* SQ_ERI_G_S_F_S_b = scr.getNewMemPos(150);
  Double* SQ_ERI_F_S_H_S_b = scr.getNewMemPos(210);
  Double* SQ_ERI_F_S_G_S_b = scr.getNewMemPos(150);
  Double* SQ_ERI_F_S_F_S_b = scr.getNewMemPos(100);
  Double* SQ_ERI_H_S_K_S_cc = scr.getNewMemPos(756);
  Double* SQ_ERI_H_S_I_S_cc = scr.getNewMemPos(588);
  Double* SQ_ERI_H_S_H_S_cc = scr.getNewMemPos(441);
  Double* SQ_ERI_G_S_K_S_cc = scr.getNewMemPos(540);
  Double* SQ_ERI_G_S_I_S_cc = scr.getNewMemPos(420);
  Double* SQ_ERI_G_S_H_S_cc = scr.getNewMemPos(315);
  Double* SQ_ERI_F_S_K_S_cc = scr.getNewMemPos(360);
  Double* SQ_ERI_F_S_I_S_cc = scr.getNewMemPos(280);
  Double* SQ_ERI_F_S_H_S_cc = scr.getNewMemPos(210);
  Double* SQ_ERI_H_S_H_S_c = scr.getNewMemPos(441);
  Double* SQ_ERI_H_S_G_S_c = scr.getNewMemPos(315);
  Double* SQ_ERI_H_S_F_S_c = scr.getNewMemPos(210);
  Double* SQ_ERI_G_S_F_S_c = scr.getNewMemPos(150);
  Double* SQ_ERI_F_S_F_S_c = scr.getNewMemPos(100);
  Double* SQ_ERI_H_S_F_S = scr.getNewMemPos(210);
  Double* SQ_ERI_H_S_D_S = scr.getNewMemPos(126);
  Double* SQ_ERI_H_S_P_S = scr.getNewMemPos(63);
  Double* SQ_ERI_G_S_P_S = scr.getNewMemPos(45);
  Double* SQ_ERI_F_S_P_S = scr.getNewMemPos(30);
  Double* SQ_ERI_K_S_H_S_ab = scr.getNewMemPos(756);
  Double* SQ_ERI_K_S_G_S_ab = scr.getNewMemPos(540);
  Double* SQ_ERI_K_S_F_S_ab = scr.getNewMemPos(360);
  Double* SQ_ERI_I_S_H_S_ab = scr.getNewMemPos(588);
  Double* SQ_ERI_I_S_G_S_ab = scr.getNewMemPos(420);
  Double* SQ_ERI_I_S_F_S_ab = scr.getNewMemPos(280);
  Double* SQ_ERI_H_S_H_S_ab = scr.getNewMemPos(441);
  Double* SQ_ERI_H_S_G_S_ab = scr.getNewMemPos(315);
  Double* SQ_ERI_H_S_F_S_ab = scr.getNewMemPos(210);
  Double* SQ_ERI_G_S_H_S_ab = scr.getNewMemPos(315);
  Double* SQ_ERI_G_S_G_S_ab = scr.getNewMemPos(225);
  Double* SQ_ERI_G_S_F_S_ab = scr.getNewMemPos(150);
  Double* SQ_ERI_D_S_H_S_b = scr.getNewMemPos(126);
  Double* SQ_ERI_D_S_G_S_b = scr.getNewMemPos(90);
  Double* SQ_ERI_D_S_F_S_b = scr.getNewMemPos(60);
  Double* SQ_ERI_I_S_I_S_bc = scr.getNewMemPos(784);
  Double* SQ_ERI_I_S_H_S_bc = scr.getNewMemPos(588);
  Double* SQ_ERI_I_S_G_S_bc = scr.getNewMemPos(420);
  Double* SQ_ERI_H_S_I_S_bc = scr.getNewMemPos(588);
  Double* SQ_ERI_H_S_H_S_bc = scr.getNewMemPos(441);
  Double* SQ_ERI_H_S_G_S_bc = scr.getNewMemPos(315);
  Double* SQ_ERI_G_S_I_S_bc = scr.getNewMemPos(420);
  Double* SQ_ERI_G_S_H_S_bc = scr.getNewMemPos(315);
  Double* SQ_ERI_G_S_G_S_bc = scr.getNewMemPos(225);
  Double* SQ_ERI_F_S_I_S_bc = scr.getNewMemPos(280);
  Double* SQ_ERI_F_S_H_S_bc = scr.getNewMemPos(210);
  Double* SQ_ERI_F_S_G_S_bc = scr.getNewMemPos(150);
  Double* SQ_ERI_I_S_G_S_b = scr.getNewMemPos(420);
  Double* SQ_ERI_I_S_F_S_b = scr.getNewMemPos(280);
  Double* SQ_ERI_I_S_D_S_b = scr.getNewMemPos(168);
  Double* SQ_ERI_H_S_D_S_b = scr.getNewMemPos(126);
  Double* SQ_ERI_G_S_D_S_b = scr.getNewMemPos(90);
  Double* SQ_ERI_F_S_D_S_b = scr.getNewMemPos(60);
  Double* SQ_ERI_K_S_H_S_bb = scr.getNewMemPos(756);
  Double* SQ_ERI_K_S_G_S_bb = scr.getNewMemPos(540);
  Double* SQ_ERI_K_S_F_S_bb = scr.getNewMemPos(360);
  Double* SQ_ERI_I_S_H_S_bb = scr.getNewMemPos(588);
  Double* SQ_ERI_I_S_G_S_bb = scr.getNewMemPos(420);
  Double* SQ_ERI_I_S_F_S_bb = scr.getNewMemPos(280);
  Double* SQ_ERI_H_S_H_S_bb = scr.getNewMemPos(441);
  Double* SQ_ERI_H_S_G_S_bb = scr.getNewMemPos(315);
  Double* SQ_ERI_H_S_F_S_bb = scr.getNewMemPos(210);
  Double* SQ_ERI_G_S_H_S_bb = scr.getNewMemPos(315);
  Double* SQ_ERI_G_S_G_S_bb = scr.getNewMemPos(225);
  Double* SQ_ERI_G_S_F_S_bb = scr.getNewMemPos(150);
  Double* SQ_ERI_F_S_H_S_bb = scr.getNewMemPos(210);
  Double* SQ_ERI_F_S_G_S_bb = scr.getNewMemPos(150);
  Double* SQ_ERI_F_S_F_S_bb = scr.getNewMemPos(100);

  //
  // declare VRR results in array form, for the case that VRR and 
  // contraction is split, so contraction will be done after VRR 
  //
  Double* SQ_ERI_F_S_H_S_vrr_array = scr.getNewMemPos(210);
  Double* SQ_ERI_F_S_G_S_vrr_array = scr.getNewMemPos(150);
  Double* SQ_ERI_F_S_F_S_vrr_array = scr.getNewMemPos(100);
  Double* SQ_ERI_H_S_H_S_vrr_array = scr.getNewMemPos(441);
  Double* SQ_ERI_H_S_G_S_vrr_array = scr.getNewMemPos(315);
  Double* SQ_ERI_H_S_F_S_vrr_array = scr.getNewMemPos(210);
  Double* SQ_ERI_G_S_H_S_vrr_array = scr.getNewMemPos(315);
  Double* SQ_ERI_G_S_G_S_vrr_array = scr.getNewMemPos(225);
  Double* SQ_ERI_G_S_F_S_vrr_array = scr.getNewMemPos(150);
  Double* SQ_ERI_D_S_H_S_vrr_array = scr.getNewMemPos(126);
  Double* SQ_ERI_D_S_G_S_vrr_array = scr.getNewMemPos(90);
  Double* SQ_ERI_D_S_F_S_vrr_array = scr.getNewMemPos(60);
  Double* SQ_ERI_G_S_I_S_vrr_array = scr.getNewMemPos(420);
  Double* SQ_ERI_F_S_I_S_vrr_array = scr.getNewMemPos(280);
  Double* SQ_ERI_G_S_D_S_vrr_array = scr.getNewMemPos(90);
  Double* SQ_ERI_F_S_D_S_vrr_array = scr.getNewMemPos(60);
  Double* SQ_ERI_K_S_H_S_vrr_array = scr.getNewMemPos(756);
  Double* SQ_ERI_K_S_G_S_vrr_array = scr.getNewMemPos(540);
  Double* SQ_ERI_K_S_F_S_vrr_array = scr.getNewMemPos(360);
  Double* SQ_ERI_I_S_H_S_vrr_array = scr.getNewMemPos(588);
  Double* SQ_ERI_I_S_G_S_vrr_array = scr.getNewMemPos(420);
  Double* SQ_ERI_I_S_F_S_vrr_array = scr.getNewMemPos(280);
  Double* SQ_ERI_P_S_H_S_vrr_array = scr.getNewMemPos(63);
  Double* SQ_ERI_P_S_G_S_vrr_array = scr.getNewMemPos(45);
  Double* SQ_ERI_P_S_F_S_vrr_array = scr.getNewMemPos(30);
  Double* SQ_ERI_I_S_I_S_vrr_array = scr.getNewMemPos(784);
  Double* SQ_ERI_H_S_I_S_vrr_array = scr.getNewMemPos(588);
  Double* SQ_ERI_I_S_D_S_vrr_array = scr.getNewMemPos(168);
  Double* SQ_ERI_H_S_D_S_vrr_array = scr.getNewMemPos(126);
  Double* SQ_ERI_D_S_I_S_vrr_array = scr.getNewMemPos(168);
  Double* SQ_ERI_D_S_D_S_vrr_array = scr.getNewMemPos(36);
  Double* SQ_ERI_H_S_K_S_vrr_array = scr.getNewMemPos(756);
  Double* SQ_ERI_G_S_K_S_vrr_array = scr.getNewMemPos(540);
  Double* SQ_ERI_F_S_K_S_vrr_array = scr.getNewMemPos(360);
  Double* SQ_ERI_H_S_P_S_vrr_array = scr.getNewMemPos(63);
  Double* SQ_ERI_G_S_P_S_vrr_array = scr.getNewMemPos(45);
  Double* SQ_ERI_F_S_P_S_vrr_array = scr.getNewMemPos(30);

  // initialize the significance check for VRR part 
  // this will determine that whether we skip the following part 
  bool isSignificant = false;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double onedz = iexp[ip2];
    Double zeta  = 1.0E0/onedz;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
    Double ic2   = icoe[ip2];
    Double fbra  = ifac[ip2];
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    for(UInt jp2=0; jp2<jnp2; jp2++) {
      Double onede = jexp[jp2];
      Double eta   = 1.0E0/onede;
      Double ediff = jexpdiff[jp2];
      Double gamma = 0.5E0*(eta+ediff);
      Double delta = 0.5E0*(eta-ediff);
      Double jc2   = jcoe[jp2];
      Double fket  = jfac[jp2];
      Double pref      = fbra*fket;
      Double prefactor = ic2*jc2*pref;

      // 
      // here below the code is performing significance test for integrals on
      // primitive integrals. Here we use the overlap integrals to roughly 
      // estimate the order of the result integrals
      // the threshold value should be for primitive function quartet, we compare
      // the value against machine precision for significance test
      // 
      Double I_ERI_S_S_S_S_vrr_IntegralTest = pref;
      if (fabs(ic2*jc2)>1.0E0) {
        I_ERI_S_S_S_S_vrr_IntegralTest = prefactor;
      }

      // test the integrals with the pMax, which is the maximum value
      // of the corresponding density matrix block(or it may be maximum
      // value pair of the corresponding density matrix block)
      if(fabs(I_ERI_S_S_S_S_vrr_IntegralTest*pMax)<THRESHOLD_MATH) continue;
      isSignificant = true;


      UInt offsetQ  = 3*jp2;
      Double QX    = Q[offsetQ  ];
      Double QY    = Q[offsetQ+1];
      Double QZ    = Q[offsetQ+2];
      Double rho   = 1.0E0/(onedz+onede);
      Double sqrho = sqrt(rho);
      Double PQ2   = (PX-QX)*(PX-QX)+(PY-QY)*(PY-QY)+(PZ-QZ)*(PZ-QZ);
      Double u     = rho*PQ2;
      if (withErfR12) u = PQ2/(1.0E0/(omega*omega)+1.0E0/rho);
      Double squ   = sqrt(u);
      Double QCX   = QX - C[0];
      Double QCY   = QY - C[1];
      Double QCZ   = QZ - C[2];
      Double WX    = rho*(PX*onede + QX*onedz);
      Double WY    = rho*(PY*onede + QY*onedz);
      Double WZ    = rho*(PZ*onede + QZ*onedz);
      Double oned2k= 0.5E0*rho*onede*onedz;
      Double WPX   = WX - PX;
      Double WPY   = WY - PY;
      Double WPZ   = WZ - PZ;
      Double rhod2zsq = rho*oned2z*onedz;
      Double WQX   = WX - QX;
      Double WQY   = WY - QY;
      Double WQZ   = WZ - QZ;
      Double oned2e= 0.5E0*onede;
      Double rhod2esq= rho*oned2e*onede;


      //
      //
      // now here for maxM>0 to compute the infamous incomplete Gamma function f_{m}(u)
      // the implementation is divided in two situations:
      // 1  if u <=1.8; use power series to get f_{Mmax}(u), then use down recursive
      //    relation to get the rest of incomplete Gamma functions;
      // 2  for u >1.8 and M <= 10 we calculate erf(u), then use up recursive
      //    relation to calculate the rest of results
      // 3  for u> 1.8 and M >  10 we calculate f_{Mmax}(u) then use down 
      //    recursive relation to get rest of incomplete Gamma functions 
      // The above procedure is tested for u between 0 to 40 with step length 1.0E-6
      // (or 1.0E-5 for float double data), for up recursive relation it shows the error
      // within 1.0E-12 (for M_limit = 12 or error within 1.0E-6 for float type of data
      // For the polynomial expansion and down recursive procedure the error is within 
      // 1.0E-14. All of the testing details please refer to the fmt_test folder
      // 
      // There's one thing need to note for up recursive process. We found that the up
      // recursive procedure is only stable for maxM<=10 and u>1.8 with double
      // precision data, single precision data will lose accuracy quickly so the result
      // for single precision calculation is not doable. Therefore if the "WITH_SINGLE_PRECISION"
      // is defined, then for erf function calculation as well as up recursive
      // process we will use the double type of data
      // 
      //

      Double I_ERI_S_S_S_S_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M1_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M2_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M3_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M4_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M5_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M6_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M7_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M8_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M9_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M10_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M11_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M12_vrr  = 0.0E0;

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER59;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER57*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER55*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER53*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER51*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER49*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER47*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER45*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER43*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER41*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER39*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER37*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER35*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER33*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER31*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER29*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = 1.0E0+u2*ONEOVER27*I_ERI_S_S_S_S_M12_vrr;
        I_ERI_S_S_S_S_M12_vrr = ONEOVER25*I_ERI_S_S_S_S_M12_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M12_vrr  = f*I_ERI_S_S_S_S_M12_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ERI_S_S_S_S_M11_vrr  = ONEOVER23*(u2*I_ERI_S_S_S_S_M12_vrr+f);
        I_ERI_S_S_S_S_M10_vrr  = ONEOVER21*(u2*I_ERI_S_S_S_S_M11_vrr+f);
        I_ERI_S_S_S_S_M9_vrr  = ONEOVER19*(u2*I_ERI_S_S_S_S_M10_vrr+f);
        I_ERI_S_S_S_S_M8_vrr  = ONEOVER17*(u2*I_ERI_S_S_S_S_M9_vrr+f);
        I_ERI_S_S_S_S_M7_vrr  = ONEOVER15*(u2*I_ERI_S_S_S_S_M8_vrr+f);
        I_ERI_S_S_S_S_M6_vrr  = ONEOVER13*(u2*I_ERI_S_S_S_S_M7_vrr+f);
        I_ERI_S_S_S_S_M5_vrr  = ONEOVER11*(u2*I_ERI_S_S_S_S_M6_vrr+f);
        I_ERI_S_S_S_S_M4_vrr  = ONEOVER9*(u2*I_ERI_S_S_S_S_M5_vrr+f);
        I_ERI_S_S_S_S_M3_vrr  = ONEOVER7*(u2*I_ERI_S_S_S_S_M4_vrr+f);
        I_ERI_S_S_S_S_M2_vrr  = ONEOVER5*(u2*I_ERI_S_S_S_S_M3_vrr+f);
        I_ERI_S_S_S_S_M1_vrr  = ONEOVER3*(u2*I_ERI_S_S_S_S_M2_vrr+f);
        I_ERI_S_S_S_S_vrr  = ONEOVER1*(u2*I_ERI_S_S_S_S_M1_vrr+f);

      }else{

        //
        // now here for maxM>M_limit
        // use external function to calculate f_{Mmax}(t)
        // then use down recursive relation to get others
        //

        // calculate (SS|SS)^{Mmax} with incomplete gamma function
        // currently we use boost library for calculation
        if (fabs(u)<THRESHOLD_MATH) {
          I_ERI_S_S_S_S_M12_vrr = 1.0E0/(2.0E0*12+1.0E0);
        }else{
          I_ERI_S_S_S_S_M12_vrr = (0.5E0/squ)*boost::math::tgamma_lower(12+0.5E0,u);
          Double oneOveru = 1.0E0/u;
          for(UInt i=0; i<12; i++) {
            I_ERI_S_S_S_S_M12_vrr = I_ERI_S_S_S_S_M12_vrr*oneOveru;
          }
        }
        Double f = TWOOVERSQRTPI*prefactor*sqrho;
        I_ERI_S_S_S_S_M12_vrr *= f;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        Double u2 = 2.0E0*u;
        Double eu = exp(-u);
        f = f*eu;
        I_ERI_S_S_S_S_M11_vrr  = ONEOVER23*(u2*I_ERI_S_S_S_S_M12_vrr+f);
        I_ERI_S_S_S_S_M10_vrr  = ONEOVER21*(u2*I_ERI_S_S_S_S_M11_vrr+f);
        I_ERI_S_S_S_S_M9_vrr  = ONEOVER19*(u2*I_ERI_S_S_S_S_M10_vrr+f);
        I_ERI_S_S_S_S_M8_vrr  = ONEOVER17*(u2*I_ERI_S_S_S_S_M9_vrr+f);
        I_ERI_S_S_S_S_M7_vrr  = ONEOVER15*(u2*I_ERI_S_S_S_S_M8_vrr+f);
        I_ERI_S_S_S_S_M6_vrr  = ONEOVER13*(u2*I_ERI_S_S_S_S_M7_vrr+f);
        I_ERI_S_S_S_S_M5_vrr  = ONEOVER11*(u2*I_ERI_S_S_S_S_M6_vrr+f);
        I_ERI_S_S_S_S_M4_vrr  = ONEOVER9*(u2*I_ERI_S_S_S_S_M5_vrr+f);
        I_ERI_S_S_S_S_M3_vrr  = ONEOVER7*(u2*I_ERI_S_S_S_S_M4_vrr+f);
        I_ERI_S_S_S_S_M2_vrr  = ONEOVER5*(u2*I_ERI_S_S_S_S_M3_vrr+f);
        I_ERI_S_S_S_S_M1_vrr  = ONEOVER3*(u2*I_ERI_S_S_S_S_M2_vrr+f);
        I_ERI_S_S_S_S_vrr  = ONEOVER1*(u2*I_ERI_S_S_S_S_M1_vrr+f);
      }


      // now scale the bottom integral if oper in erf(r12)/r12 form
      if (withErfR12) {
        Double erfPref0   = 1.0E0+rho/(omega*omega);
        Double erfPref1   = 1.0E0/erfPref0;
        Double erfp       = sqrt(erfPref1);
        Double erfp2      = erfp*erfp;
        Double erfPref_1  = erfp;
        I_ERI_S_S_S_S_vrr = I_ERI_S_S_S_S_vrr*erfPref_1;
        Double erfPref_3 = erfPref_1*erfp2;
        Double erfPref_5 = erfPref_3*erfp2;
        Double erfPref_7 = erfPref_5*erfp2;
        Double erfPref_9 = erfPref_7*erfp2;
        Double erfPref_11 = erfPref_9*erfp2;
        Double erfPref_13 = erfPref_11*erfp2;
        Double erfPref_15 = erfPref_13*erfp2;
        Double erfPref_17 = erfPref_15*erfp2;
        Double erfPref_19 = erfPref_17*erfp2;
        Double erfPref_21 = erfPref_19*erfp2;
        Double erfPref_23 = erfPref_21*erfp2;
        Double erfPref_25 = erfPref_23*erfp2;
        I_ERI_S_S_S_S_M1_vrr = I_ERI_S_S_S_S_M1_vrr*erfPref_3;
        I_ERI_S_S_S_S_M2_vrr = I_ERI_S_S_S_S_M2_vrr*erfPref_5;
        I_ERI_S_S_S_S_M3_vrr = I_ERI_S_S_S_S_M3_vrr*erfPref_7;
        I_ERI_S_S_S_S_M4_vrr = I_ERI_S_S_S_S_M4_vrr*erfPref_9;
        I_ERI_S_S_S_S_M5_vrr = I_ERI_S_S_S_S_M5_vrr*erfPref_11;
        I_ERI_S_S_S_S_M6_vrr = I_ERI_S_S_S_S_M6_vrr*erfPref_13;
        I_ERI_S_S_S_S_M7_vrr = I_ERI_S_S_S_S_M7_vrr*erfPref_15;
        I_ERI_S_S_S_S_M8_vrr = I_ERI_S_S_S_S_M8_vrr*erfPref_17;
        I_ERI_S_S_S_S_M9_vrr = I_ERI_S_S_S_S_M9_vrr*erfPref_19;
        I_ERI_S_S_S_S_M10_vrr = I_ERI_S_S_S_S_M10_vrr*erfPref_21;
        I_ERI_S_S_S_S_M11_vrr = I_ERI_S_S_S_S_M11_vrr*erfPref_23;
        I_ERI_S_S_S_S_M12_vrr = I_ERI_S_S_S_S_M12_vrr*erfPref_25;
      }

      hgp_os_eri_f_d_f_d_d2_vrr_1(oned2k,oned2z,rhod2zsq,WPX,WPY,WPZ,oned2e,rhod2esq,WQX,WQY,WQZ,PAX,PAY,PAZ,QCX,QCY,QCZ,alpha,beta,gamma,delta,I_ERI_S_S_S_S_vrr,I_ERI_S_S_S_S_M1_vrr,I_ERI_S_S_S_S_M2_vrr,I_ERI_S_S_S_S_M3_vrr,I_ERI_S_S_S_S_M4_vrr,I_ERI_S_S_S_S_M5_vrr,I_ERI_S_S_S_S_M6_vrr,I_ERI_S_S_S_S_M7_vrr,I_ERI_S_S_S_S_M8_vrr,I_ERI_S_S_S_S_M9_vrr,I_ERI_S_S_S_S_M10_vrr,I_ERI_S_S_S_S_M11_vrr,I_ERI_S_S_S_S_M12_vrr,&SQ_ERI_D_S_D_S_vrr_array[0],&SQ_ERI_P_S_F_S_vrr_array[0],&SQ_ERI_F_S_P_S_vrr_array[0],&SQ_ERI_D_S_F_S_vrr_array[0],&SQ_ERI_F_S_D_S_vrr_array[0],&SQ_ERI_P_S_G_S_vrr_array[0],&SQ_ERI_G_S_P_S_vrr_array[0],&SQ_ERI_F_S_F_S_vrr_array[0],&SQ_ERI_D_S_G_S_vrr_array[0],&SQ_ERI_G_S_D_S_vrr_array[0],&SQ_ERI_P_S_H_S_vrr_array[0],&SQ_ERI_H_S_P_S_vrr_array[0],&SQ_ERI_F_S_G_S_vrr_array[0],&SQ_ERI_G_S_F_S_vrr_array[0],&SQ_ERI_D_S_H_S_vrr_array[0],&SQ_ERI_H_S_D_S_vrr_array[0],&SQ_ERI_G_S_G_S_vrr_array[0],&SQ_ERI_F_S_H_S_vrr_array[0],&SQ_ERI_H_S_F_S_vrr_array[0],&SQ_ERI_D_S_I_S_vrr_array[0],&SQ_ERI_I_S_D_S_vrr_array[0],&SQ_ERI_G_S_H_S_vrr_array[0],&SQ_ERI_H_S_G_S_vrr_array[0],&SQ_ERI_F_S_I_S_vrr_array[0],&SQ_ERI_I_S_F_S_vrr_array[0],&SQ_ERI_H_S_H_S_vrr_array[0],&SQ_ERI_G_S_I_S_vrr_array[0],&SQ_ERI_I_S_G_S_vrr_array[0],&SQ_ERI_F_S_K_S_vrr_array[0],&SQ_ERI_K_S_F_S_vrr_array[0],&SQ_ERI_H_S_I_S_vrr_array[0],&SQ_ERI_I_S_H_S_vrr_array[0],&SQ_ERI_G_S_K_S_vrr_array[0],&SQ_ERI_K_S_G_S_vrr_array[0],&SQ_ERI_I_S_I_S_vrr_array[0],&SQ_ERI_H_S_K_S_vrr_array[0],&SQ_ERI_K_S_H_S_vrr_array[0]);



      hgp_os_eri_f_d_f_d_d2_vrrcont_1(alpha,beta,gamma,delta,&SQ_ERI_F_S_H_S_vrr_array[0],&SQ_ERI_F_S_G_S_vrr_array[0],&SQ_ERI_F_S_F_S_vrr_array[0],&SQ_ERI_H_S_H_S_vrr_array[0],&SQ_ERI_H_S_G_S_vrr_array[0],&SQ_ERI_H_S_F_S_vrr_array[0],&SQ_ERI_G_S_H_S_vrr_array[0],&SQ_ERI_G_S_G_S_vrr_array[0],&SQ_ERI_F_S_H_S[0],&SQ_ERI_F_S_G_S[0],&SQ_ERI_F_S_F_S[0],&SQ_ERI_H_S_H_S_a[0],&SQ_ERI_H_S_G_S_a[0],&SQ_ERI_H_S_F_S_a[0],&SQ_ERI_G_S_H_S_a[0],&SQ_ERI_G_S_G_S_a[0],&SQ_ERI_G_S_H_S_c[0],&SQ_ERI_G_S_G_S_c[0],&SQ_ERI_F_S_H_S_c[0],&SQ_ERI_F_S_G_S_c[0],&SQ_ERI_G_S_G_S[0],&SQ_ERI_H_S_H_S_aa[0],&SQ_ERI_H_S_G_S_aa[0],&SQ_ERI_H_S_F_S_aa[0],&SQ_ERI_F_S_H_S_a[0],&SQ_ERI_F_S_G_S_a[0],&SQ_ERI_F_S_F_S_a[0],&SQ_ERI_H_S_H_S_ac[0],&SQ_ERI_H_S_G_S_ac[0],&SQ_ERI_G_S_H_S_ac[0],&SQ_ERI_G_S_G_S_ac[0],&SQ_ERI_H_S_H_S_b[0],&SQ_ERI_H_S_G_S_b[0],&SQ_ERI_H_S_F_S_b[0],&SQ_ERI_G_S_H_S_b[0],&SQ_ERI_G_S_G_S_b[0],&SQ_ERI_F_S_H_S_b[0],&SQ_ERI_F_S_G_S_b[0],&SQ_ERI_F_S_F_S_b[0],&SQ_ERI_H_S_H_S_cc[0],&SQ_ERI_G_S_H_S_cc[0],&SQ_ERI_F_S_H_S_cc[0],&SQ_ERI_H_S_H_S_c[0],&SQ_ERI_H_S_G_S_c[0],&SQ_ERI_H_S_F_S_c[0],&SQ_ERI_F_S_F_S_c[0],&SQ_ERI_H_S_F_S[0],&SQ_ERI_H_S_H_S_ab[0],&SQ_ERI_H_S_G_S_ab[0],&SQ_ERI_H_S_F_S_ab[0],&SQ_ERI_G_S_H_S_ab[0],&SQ_ERI_G_S_G_S_ab[0],&SQ_ERI_H_S_H_S_bc[0],&SQ_ERI_H_S_G_S_bc[0],&SQ_ERI_G_S_H_S_bc[0],&SQ_ERI_G_S_G_S_bc[0],&SQ_ERI_F_S_H_S_bc[0],&SQ_ERI_F_S_G_S_bc[0],&SQ_ERI_H_S_H_S_bb[0],&SQ_ERI_H_S_G_S_bb[0],&SQ_ERI_H_S_F_S_bb[0],&SQ_ERI_G_S_H_S_bb[0],&SQ_ERI_G_S_G_S_bb[0],&SQ_ERI_F_S_H_S_bb[0],&SQ_ERI_F_S_G_S_bb[0],&SQ_ERI_F_S_F_S_bb[0]);

      hgp_os_eri_f_d_f_d_d2_vrrcont_2(alpha,beta,gamma,delta,&SQ_ERI_G_S_F_S_vrr_array[0],&SQ_ERI_D_S_H_S_vrr_array[0],&SQ_ERI_D_S_G_S_vrr_array[0],&SQ_ERI_D_S_F_S_vrr_array[0],&SQ_ERI_G_S_I_S_vrr_array[0],&SQ_ERI_F_S_I_S_vrr_array[0],&SQ_ERI_G_S_D_S_vrr_array[0],&SQ_ERI_F_S_D_S_vrr_array[0],&SQ_ERI_K_S_H_S_vrr_array[0],&SQ_ERI_K_S_G_S_vrr_array[0],&SQ_ERI_K_S_F_S_vrr_array[0],&SQ_ERI_I_S_H_S_vrr_array[0],&SQ_ERI_I_S_G_S_vrr_array[0],&SQ_ERI_I_S_F_S_vrr_array[0],&SQ_ERI_G_S_F_S_a[0],&SQ_ERI_D_S_H_S[0],&SQ_ERI_D_S_G_S[0],&SQ_ERI_D_S_F_S[0],&SQ_ERI_G_S_I_S_c[0],&SQ_ERI_F_S_I_S_c[0],&SQ_ERI_G_S_F_S[0],&SQ_ERI_G_S_D_S[0],&SQ_ERI_F_S_D_S[0],&SQ_ERI_K_S_H_S_aa[0],&SQ_ERI_K_S_G_S_aa[0],&SQ_ERI_K_S_F_S_aa[0],&SQ_ERI_I_S_H_S_aa[0],&SQ_ERI_I_S_G_S_aa[0],&SQ_ERI_I_S_F_S_aa[0],&SQ_ERI_I_S_H_S_ac[0],&SQ_ERI_I_S_G_S_ac[0],&SQ_ERI_G_S_I_S_ac[0],&SQ_ERI_I_S_G_S_a[0],&SQ_ERI_I_S_F_S_a[0],&SQ_ERI_G_S_D_S_a[0],&SQ_ERI_D_S_H_S_c[0],&SQ_ERI_D_S_G_S_c[0],&SQ_ERI_G_S_F_S_b[0],&SQ_ERI_G_S_I_S_cc[0],&SQ_ERI_F_S_I_S_cc[0],&SQ_ERI_G_S_F_S_c[0],&SQ_ERI_K_S_H_S_ab[0],&SQ_ERI_K_S_G_S_ab[0],&SQ_ERI_K_S_F_S_ab[0],&SQ_ERI_I_S_H_S_ab[0],&SQ_ERI_I_S_G_S_ab[0],&SQ_ERI_I_S_F_S_ab[0],&SQ_ERI_G_S_F_S_ab[0],&SQ_ERI_D_S_H_S_b[0],&SQ_ERI_D_S_G_S_b[0],&SQ_ERI_D_S_F_S_b[0],&SQ_ERI_I_S_H_S_bc[0],&SQ_ERI_I_S_G_S_bc[0],&SQ_ERI_G_S_I_S_bc[0],&SQ_ERI_F_S_I_S_bc[0],&SQ_ERI_I_S_G_S_b[0],&SQ_ERI_I_S_F_S_b[0],&SQ_ERI_G_S_D_S_b[0],&SQ_ERI_F_S_D_S_b[0],&SQ_ERI_K_S_H_S_bb[0],&SQ_ERI_K_S_G_S_bb[0],&SQ_ERI_K_S_F_S_bb[0],&SQ_ERI_I_S_H_S_bb[0],&SQ_ERI_I_S_G_S_bb[0],&SQ_ERI_I_S_F_S_bb[0],&SQ_ERI_G_S_F_S_bb[0]);

      hgp_os_eri_f_d_f_d_d2_vrrcont_3(alpha,beta,gamma,delta,&SQ_ERI_P_S_H_S_vrr_array[0],&SQ_ERI_P_S_G_S_vrr_array[0],&SQ_ERI_P_S_F_S_vrr_array[0],&SQ_ERI_I_S_I_S_vrr_array[0],&SQ_ERI_H_S_I_S_vrr_array[0],&SQ_ERI_I_S_D_S_vrr_array[0],&SQ_ERI_H_S_D_S_vrr_array[0],&SQ_ERI_D_S_I_S_vrr_array[0],&SQ_ERI_D_S_D_S_vrr_array[0],&SQ_ERI_H_S_K_S_vrr_array[0],&SQ_ERI_G_S_K_S_vrr_array[0],&SQ_ERI_F_S_K_S_vrr_array[0],&SQ_ERI_H_S_P_S_vrr_array[0],&SQ_ERI_G_S_P_S_vrr_array[0],&SQ_ERI_F_S_P_S_vrr_array[0],&SQ_ERI_P_S_H_S[0],&SQ_ERI_P_S_G_S[0],&SQ_ERI_P_S_F_S[0],&SQ_ERI_I_S_I_S_ac[0],&SQ_ERI_H_S_I_S_ac[0],&SQ_ERI_I_S_D_S_a[0],&SQ_ERI_H_S_D_S_a[0],&SQ_ERI_D_S_I_S_c[0],&SQ_ERI_D_S_D_S[0],&SQ_ERI_H_S_K_S_cc[0],&SQ_ERI_H_S_I_S_cc[0],&SQ_ERI_G_S_K_S_cc[0],&SQ_ERI_F_S_K_S_cc[0],&SQ_ERI_H_S_D_S[0],&SQ_ERI_H_S_P_S[0],&SQ_ERI_G_S_P_S[0],&SQ_ERI_F_S_P_S[0],&SQ_ERI_I_S_I_S_bc[0],&SQ_ERI_H_S_I_S_bc[0],&SQ_ERI_I_S_D_S_b[0],&SQ_ERI_H_S_D_S_b[0]);


    }
  }

  /************************************************************
   * let's see the significance test result. if VRR result is
   * insignificant, there's no need to do following codes
   ************************************************************/
  if (! isSignificant) return;

  /************************************************************
   * declare the HRR1 result shell quartets in array form
   ************************************************************/
  Double* SQ_ERI_F_S_F_D = scr.getNewMemPos(600);
  Double* SQ_ERI_H_S_F_D_a = scr.getNewMemPos(1260);
  Double* SQ_ERI_G_S_F_D_a = scr.getNewMemPos(900);
  Double* SQ_ERI_D_S_F_D = scr.getNewMemPos(360);
  Double* SQ_ERI_G_S_G_D_c = scr.getNewMemPos(1350);
  Double* SQ_ERI_F_S_G_D_c = scr.getNewMemPos(900);
  Double* SQ_ERI_G_S_D_D = scr.getNewMemPos(540);
  Double* SQ_ERI_F_S_D_D = scr.getNewMemPos(360);
  Double* SQ_ERI_K_S_F_D_aa = scr.getNewMemPos(2160);
  Double* SQ_ERI_I_S_F_D_aa = scr.getNewMemPos(1680);
  Double* SQ_ERI_H_S_F_D_aa = scr.getNewMemPos(1260);
  Double* SQ_ERI_F_S_F_D_a = scr.getNewMemPos(600);
  Double* SQ_ERI_P_S_F_D = scr.getNewMemPos(180);
  Double* SQ_ERI_I_S_G_D_ac = scr.getNewMemPos(2520);
  Double* SQ_ERI_H_S_G_D_ac = scr.getNewMemPos(1890);
  Double* SQ_ERI_G_S_G_D_ac = scr.getNewMemPos(1350);
  Double* SQ_ERI_I_S_D_D_a = scr.getNewMemPos(1008);
  Double* SQ_ERI_H_S_D_D_a = scr.getNewMemPos(756);
  Double* SQ_ERI_G_S_D_D_a = scr.getNewMemPos(540);
  Double* SQ_ERI_D_S_G_D_c = scr.getNewMemPos(540);
  Double* SQ_ERI_D_S_D_D = scr.getNewMemPos(216);
  Double* SQ_ERI_H_S_F_D_b = scr.getNewMemPos(1260);
  Double* SQ_ERI_G_S_F_D_b = scr.getNewMemPos(900);
  Double* SQ_ERI_F_S_F_D_b = scr.getNewMemPos(600);
  Double* SQ_ERI_H_S_H_D_cc = scr.getNewMemPos(2646);
  Double* SQ_ERI_G_S_H_D_cc = scr.getNewMemPos(1890);
  Double* SQ_ERI_F_S_H_D_cc = scr.getNewMemPos(1260);
  Double* SQ_ERI_H_S_F_D_c = scr.getNewMemPos(1260);
  Double* SQ_ERI_G_S_F_D_c = scr.getNewMemPos(900);
  Double* SQ_ERI_F_S_F_D_c = scr.getNewMemPos(600);
  Double* SQ_ERI_H_S_P_D = scr.getNewMemPos(378);
  Double* SQ_ERI_G_S_P_D = scr.getNewMemPos(270);
  Double* SQ_ERI_F_S_P_D = scr.getNewMemPos(180);
  Double* SQ_ERI_K_S_F_D_ab = scr.getNewMemPos(2160);
  Double* SQ_ERI_I_S_F_D_ab = scr.getNewMemPos(1680);
  Double* SQ_ERI_H_S_F_D_ab = scr.getNewMemPos(1260);
  Double* SQ_ERI_G_S_F_D_ab = scr.getNewMemPos(900);
  Double* SQ_ERI_D_S_F_D_b = scr.getNewMemPos(360);
  Double* SQ_ERI_I_S_G_D_bc = scr.getNewMemPos(2520);
  Double* SQ_ERI_H_S_G_D_bc = scr.getNewMemPos(1890);
  Double* SQ_ERI_G_S_G_D_bc = scr.getNewMemPos(1350);
  Double* SQ_ERI_F_S_G_D_bc = scr.getNewMemPos(900);
  Double* SQ_ERI_I_S_D_D_b = scr.getNewMemPos(1008);
  Double* SQ_ERI_H_S_D_D_b = scr.getNewMemPos(756);
  Double* SQ_ERI_G_S_D_D_b = scr.getNewMemPos(540);
  Double* SQ_ERI_F_S_D_D_b = scr.getNewMemPos(360);
  Double* SQ_ERI_K_S_F_D_bb = scr.getNewMemPos(2160);
  Double* SQ_ERI_I_S_F_D_bb = scr.getNewMemPos(1680);
  Double* SQ_ERI_H_S_F_D_bb = scr.getNewMemPos(1260);
  Double* SQ_ERI_G_S_F_D_bb = scr.getNewMemPos(900);
  Double* SQ_ERI_F_S_F_D_bb = scr.getNewMemPos(600);
  Double* SQ_ERI_G_S_F_P_ab = scr.getNewMemPos(450);
  Double* SQ_ERI_H_S_F_P_ab = scr.getNewMemPos(630);


  hgp_os_eri_f_d_f_d_d2_hrr1_1(C,D,&SQ_ERI_F_S_D_S[0],&SQ_ERI_F_S_P_S[0],&SQ_ERI_G_S_D_S[0],&SQ_ERI_G_S_P_S[0],&SQ_ERI_H_S_D_S[0],&SQ_ERI_H_S_P_S[0],&SQ_ERI_D_S_F_S[0],&SQ_ERI_D_S_D_S[0],&SQ_ERI_F_S_F_S[0],&SQ_ERI_G_S_F_S[0],&SQ_ERI_H_S_F_S[0],&SQ_ERI_P_S_G_S[0],&SQ_ERI_P_S_F_S[0],&SQ_ERI_D_S_G_S[0],&SQ_ERI_F_S_G_S[0],&SQ_ERI_G_S_G_S[0],&SQ_ERI_P_S_H_S[0],&SQ_ERI_D_S_H_S[0],&SQ_ERI_F_S_H_S[0],&SQ_ERI_G_S_F_S_a[0],&SQ_ERI_G_S_D_S_a[0],&SQ_ERI_H_S_F_S_a[0],&SQ_ERI_H_S_D_S_a[0],&SQ_ERI_I_S_F_S_a[0],&SQ_ERI_I_S_D_S_a[0],&SQ_ERI_F_S_G_S_a[0],&SQ_ERI_F_S_F_S_a[0],&SQ_ERI_G_S_G_S_a[0],&SQ_ERI_H_S_G_S_a[0],&SQ_ERI_I_S_G_S_a[0],&SQ_ERI_F_S_H_S_a[0],&SQ_ERI_G_S_H_S_a[0],&SQ_ERI_H_S_H_S_a[0],&SQ_ERI_F_S_F_S_b[0],&SQ_ERI_F_S_D_S_b[0],&SQ_ERI_G_S_F_S_b[0],&SQ_ERI_G_S_D_S_b[0],&SQ_ERI_H_S_F_S_b[0],&SQ_ERI_H_S_D_S_b[0],&SQ_ERI_I_S_F_S_b[0],&SQ_ERI_I_S_D_S_b[0],&SQ_ERI_D_S_G_S_b[0],&SQ_ERI_D_S_F_S_b[0],&SQ_ERI_F_S_G_S_b[0],&SQ_ERI_G_S_G_S_b[0],&SQ_ERI_H_S_G_S_b[0],&SQ_ERI_I_S_G_S_b[0],&SQ_ERI_D_S_H_S_b[0],&SQ_ERI_F_S_H_S_b[0],&SQ_ERI_G_S_H_S_b[0],&SQ_ERI_H_S_H_S_b[0],&SQ_ERI_F_S_G_S_c[0],&SQ_ERI_F_S_F_S_c[0],&SQ_ERI_G_S_G_S_c[0],&SQ_ERI_G_S_F_S_c[0],&SQ_ERI_H_S_G_S_c[0],&SQ_ERI_H_S_F_S_c[0],&SQ_ERI_D_S_H_S_c[0],&SQ_ERI_D_S_G_S_c[0],&SQ_ERI_F_S_H_S_c[0],&SQ_ERI_G_S_H_S_c[0],&SQ_ERI_H_S_H_S_c[0],&SQ_ERI_D_S_I_S_c[0],&SQ_ERI_F_S_I_S_c[0],&SQ_ERI_G_S_I_S_c[0],&SQ_ERI_H_S_G_S_aa[0],&SQ_ERI_H_S_F_S_aa[0],&SQ_ERI_I_S_G_S_aa[0],&SQ_ERI_I_S_F_S_aa[0],&SQ_ERI_K_S_G_S_aa[0],&SQ_ERI_K_S_F_S_aa[0],&SQ_ERI_H_S_H_S_aa[0],&SQ_ERI_I_S_H_S_aa[0],&SQ_ERI_K_S_H_S_aa[0],&SQ_ERI_G_S_G_S_ab[0],&SQ_ERI_G_S_F_S_ab[0],&SQ_ERI_H_S_G_S_ab[0],&SQ_ERI_H_S_F_S_ab[0],&SQ_ERI_F_S_P_D[0],&SQ_ERI_G_S_P_D[0],&SQ_ERI_H_S_P_D[0],&SQ_ERI_D_S_D_D[0],&SQ_ERI_F_S_D_D[0],&SQ_ERI_G_S_D_D[0],&SQ_ERI_P_S_F_D[0],&SQ_ERI_D_S_F_D[0],&SQ_ERI_F_S_F_D[0],&SQ_ERI_G_S_D_D_a[0],&SQ_ERI_H_S_D_D_a[0],&SQ_ERI_I_S_D_D_a[0],&SQ_ERI_F_S_F_D_a[0],&SQ_ERI_G_S_F_D_a[0],&SQ_ERI_H_S_F_D_a[0],&SQ_ERI_F_S_D_D_b[0],&SQ_ERI_G_S_D_D_b[0],&SQ_ERI_H_S_D_D_b[0],&SQ_ERI_I_S_D_D_b[0],&SQ_ERI_D_S_F_D_b[0],&SQ_ERI_F_S_F_D_b[0],&SQ_ERI_G_S_F_D_b[0],&SQ_ERI_H_S_F_D_b[0],&SQ_ERI_F_S_F_D_c[0],&SQ_ERI_G_S_F_D_c[0],&SQ_ERI_H_S_F_D_c[0],&SQ_ERI_D_S_G_D_c[0],&SQ_ERI_F_S_G_D_c[0],&SQ_ERI_G_S_G_D_c[0],&SQ_ERI_H_S_F_D_aa[0],&SQ_ERI_I_S_F_D_aa[0],&SQ_ERI_K_S_F_D_aa[0],&SQ_ERI_G_S_F_P_ab[0],&SQ_ERI_H_S_F_P_ab[0]);

  hgp_os_eri_f_d_f_d_d2_hrr1_2(C,D,&SQ_ERI_I_S_G_S_ab[0],&SQ_ERI_I_S_F_S_ab[0],&SQ_ERI_K_S_G_S_ab[0],&SQ_ERI_K_S_F_S_ab[0],&SQ_ERI_G_S_H_S_ab[0],&SQ_ERI_G_S_G_S_ab[0],&SQ_ERI_H_S_H_S_ab[0],&SQ_ERI_H_S_G_S_ab[0],&SQ_ERI_I_S_H_S_ab[0],&SQ_ERI_K_S_H_S_ab[0],&SQ_ERI_G_S_F_P_ab[0],&SQ_ERI_H_S_F_P_ab[0],&SQ_ERI_G_S_H_S_ac[0],&SQ_ERI_G_S_G_S_ac[0],&SQ_ERI_H_S_H_S_ac[0],&SQ_ERI_H_S_G_S_ac[0],&SQ_ERI_I_S_H_S_ac[0],&SQ_ERI_I_S_G_S_ac[0],&SQ_ERI_G_S_I_S_ac[0],&SQ_ERI_H_S_I_S_ac[0],&SQ_ERI_I_S_I_S_ac[0],&SQ_ERI_F_S_G_S_bb[0],&SQ_ERI_F_S_F_S_bb[0],&SQ_ERI_G_S_G_S_bb[0],&SQ_ERI_G_S_F_S_bb[0],&SQ_ERI_H_S_G_S_bb[0],&SQ_ERI_H_S_F_S_bb[0],&SQ_ERI_I_S_G_S_bb[0],&SQ_ERI_I_S_F_S_bb[0],&SQ_ERI_K_S_G_S_bb[0],&SQ_ERI_K_S_F_S_bb[0],&SQ_ERI_F_S_H_S_bb[0],&SQ_ERI_G_S_H_S_bb[0],&SQ_ERI_H_S_H_S_bb[0],&SQ_ERI_I_S_H_S_bb[0],&SQ_ERI_K_S_H_S_bb[0],&SQ_ERI_F_S_H_S_bc[0],&SQ_ERI_F_S_G_S_bc[0],&SQ_ERI_G_S_H_S_bc[0],&SQ_ERI_G_S_G_S_bc[0],&SQ_ERI_H_S_H_S_bc[0],&SQ_ERI_H_S_G_S_bc[0],&SQ_ERI_I_S_H_S_bc[0],&SQ_ERI_I_S_G_S_bc[0],&SQ_ERI_F_S_I_S_bc[0],&SQ_ERI_G_S_I_S_bc[0],&SQ_ERI_H_S_I_S_bc[0],&SQ_ERI_I_S_I_S_bc[0],&SQ_ERI_G_S_F_D_ab[0],&SQ_ERI_H_S_F_D_ab[0],&SQ_ERI_I_S_F_D_ab[0],&SQ_ERI_K_S_F_D_ab[0],&SQ_ERI_G_S_G_D_ac[0],&SQ_ERI_H_S_G_D_ac[0],&SQ_ERI_I_S_G_D_ac[0],&SQ_ERI_F_S_F_D_bb[0],&SQ_ERI_G_S_F_D_bb[0],&SQ_ERI_H_S_F_D_bb[0],&SQ_ERI_I_S_F_D_bb[0],&SQ_ERI_K_S_F_D_bb[0],&SQ_ERI_F_S_G_D_bc[0],&SQ_ERI_G_S_G_D_bc[0],&SQ_ERI_H_S_G_D_bc[0],&SQ_ERI_I_S_G_D_bc[0]);

  hgp_os_eri_f_d_f_d_d2_hrr1_3(C,D,&SQ_ERI_F_S_I_S_cc[0],&SQ_ERI_F_S_H_S_cc[0],&SQ_ERI_G_S_I_S_cc[0],&SQ_ERI_G_S_H_S_cc[0],&SQ_ERI_H_S_I_S_cc[0],&SQ_ERI_H_S_H_S_cc[0],&SQ_ERI_F_S_K_S_cc[0],&SQ_ERI_G_S_K_S_cc[0],&SQ_ERI_H_S_K_S_cc[0],&SQ_ERI_F_S_H_D_cc[0],&SQ_ERI_G_S_H_D_cc[0],&SQ_ERI_H_S_H_D_cc[0]);



  /************************************************************
   * declare the HRR2 result shell quartets in array form
   ************************************************************/
  Double* SQ_ERI_H_D_F_D_aa = scr.getNewMemPos(7560);
  Double* SQ_ERI_F_D_F_D_a = scr.getNewMemPos(3600);
  Double* SQ_ERI_P_D_F_D = scr.getNewMemPos(1080);
  Double* SQ_ERI_G_F_F_D_ab = scr.getNewMemPos(9000);
  Double* SQ_ERI_G_P_F_D_a = scr.getNewMemPos(2700);
  Double* SQ_ERI_D_F_F_D_b = scr.getNewMemPos(3600);
  Double* SQ_ERI_D_P_F_D = scr.getNewMemPos(1080);
  Double* SQ_ERI_G_D_G_D_ac = scr.getNewMemPos(8100);
  Double* SQ_ERI_G_D_D_D_a = scr.getNewMemPos(3240);
  Double* SQ_ERI_D_D_G_D_c = scr.getNewMemPos(3240);
  Double* SQ_ERI_D_D_D_D = scr.getNewMemPos(1296);
  Double* SQ_ERI_F_G_F_D_bb = scr.getNewMemPos(9000);
  Double* SQ_ERI_F_D_F_D_b = scr.getNewMemPos(3600);
  Double* SQ_ERI_F_F_G_D_bc = scr.getNewMemPos(9000);
  Double* SQ_ERI_F_F_D_D_b = scr.getNewMemPos(3600);
  Double* SQ_ERI_F_P_G_D_c = scr.getNewMemPos(2700);
  Double* SQ_ERI_F_P_D_D = scr.getNewMemPos(1080);
  Double* SQ_ERI_F_D_H_D_cc = scr.getNewMemPos(7560);
  Double* SQ_ERI_F_D_F_D_c = scr.getNewMemPos(3600);
  Double* SQ_ERI_F_D_P_D = scr.getNewMemPos(1080);
  Double* SQ_ERI_D_P_G_D_c = scr.getNewMemPos(1620);
  Double* SQ_ERI_F_P_F_D_c = scr.getNewMemPos(1800);


  hgp_os_eri_f_d_f_d_d2_hrr2_1(A,B,&SQ_ERI_D_S_F_D[0],&SQ_ERI_P_S_F_D[0],&SQ_ERI_F_S_D_D[0],&SQ_ERI_D_S_D_D[0],&SQ_ERI_F_S_F_D[0],&SQ_ERI_G_S_P_D[0],&SQ_ERI_F_S_P_D[0],&SQ_ERI_G_S_D_D[0],&SQ_ERI_H_S_P_D[0],&SQ_ERI_G_S_F_D_a[0],&SQ_ERI_F_S_F_D_a[0],&SQ_ERI_H_S_D_D_a[0],&SQ_ERI_G_S_D_D_a[0],&SQ_ERI_H_S_F_D_a[0],&SQ_ERI_I_S_D_D_a[0],&SQ_ERI_F_S_F_D_b[0],&SQ_ERI_D_S_F_D_b[0],&SQ_ERI_G_S_D_D_b[0],&SQ_ERI_F_S_D_D_b[0],&SQ_ERI_G_S_F_D_b[0],&SQ_ERI_H_S_D_D_b[0],&SQ_ERI_H_S_F_D_b[0],&SQ_ERI_I_S_D_D_b[0],&SQ_ERI_F_S_G_D_c[0],&SQ_ERI_D_S_G_D_c[0],&SQ_ERI_G_S_F_D_c[0],&SQ_ERI_F_S_F_D_c[0],&SQ_ERI_D_P_F_D[0],&SQ_ERI_P_D_F_D[0],&SQ_ERI_F_P_D_D[0],&SQ_ERI_D_D_D_D[0],&SQ_ERI_F_D_P_D[0],&SQ_ERI_G_P_F_D_a[0],&SQ_ERI_F_D_F_D_a[0],&SQ_ERI_G_D_D_D_a[0],&SQ_ERI_F_D_F_D_b[0],&SQ_ERI_D_F_F_D_b[0],&SQ_ERI_F_F_D_D_b[0],&SQ_ERI_D_P_G_D_c[0],&SQ_ERI_F_P_F_D_c[0]);

  hgp_os_eri_f_d_f_d_d2_hrr2_2(A,B,&SQ_ERI_G_S_G_D_c[0],&SQ_ERI_F_S_G_D_c[0],&SQ_ERI_D_P_G_D_c[0],&SQ_ERI_H_S_F_D_c[0],&SQ_ERI_G_S_F_D_c[0],&SQ_ERI_F_P_F_D_c[0],&SQ_ERI_I_S_F_D_aa[0],&SQ_ERI_H_S_F_D_aa[0],&SQ_ERI_K_S_F_D_aa[0],&SQ_ERI_H_S_F_D_ab[0],&SQ_ERI_G_S_F_D_ab[0],&SQ_ERI_I_S_F_D_ab[0],&SQ_ERI_K_S_F_D_ab[0],&SQ_ERI_F_P_G_D_c[0],&SQ_ERI_D_D_G_D_c[0],&SQ_ERI_F_D_F_D_c[0],&SQ_ERI_H_D_F_D_aa[0],&SQ_ERI_G_F_F_D_ab[0]);

  hgp_os_eri_f_d_f_d_d2_hrr2_3(A,B,&SQ_ERI_H_S_G_D_ac[0],&SQ_ERI_G_S_G_D_ac[0],&SQ_ERI_I_S_G_D_ac[0],&SQ_ERI_G_S_F_D_bb[0],&SQ_ERI_F_S_F_D_bb[0],&SQ_ERI_H_S_F_D_bb[0],&SQ_ERI_I_S_F_D_bb[0],&SQ_ERI_K_S_F_D_bb[0],&SQ_ERI_G_D_G_D_ac[0],&SQ_ERI_F_G_F_D_bb[0]);

  hgp_os_eri_f_d_f_d_d2_hrr2_4(A,B,&SQ_ERI_G_S_G_D_bc[0],&SQ_ERI_F_S_G_D_bc[0],&SQ_ERI_H_S_G_D_bc[0],&SQ_ERI_I_S_G_D_bc[0],&SQ_ERI_G_S_H_D_cc[0],&SQ_ERI_F_S_H_D_cc[0],&SQ_ERI_H_S_H_D_cc[0],&SQ_ERI_F_F_G_D_bc[0],&SQ_ERI_F_D_H_D_cc[0]);



  hgp_os_eri_f_d_f_d_d2_deriv_1(&SQ_ERI_H_D_F_D_aa[0],&SQ_ERI_F_D_F_D_a[0],&SQ_ERI_P_D_F_D[0],&SQ_ERI_G_F_F_D_ab[0],&SQ_ERI_G_P_F_D_a[0],&SQ_ERI_D_F_F_D_b[0],&SQ_ERI_D_P_F_D[0],abcd);

  hgp_os_eri_f_d_f_d_d2_deriv_2(&SQ_ERI_G_F_F_D_ab[0],&SQ_ERI_G_P_F_D_a[0],&SQ_ERI_D_F_F_D_b[0],&SQ_ERI_D_P_F_D[0],&SQ_ERI_G_D_G_D_ac[0],&SQ_ERI_G_D_D_D_a[0],&SQ_ERI_D_D_G_D_c[0],&SQ_ERI_D_D_D_D[0],abcd);

  hgp_os_eri_f_d_f_d_d2_deriv_3(&SQ_ERI_F_G_F_D_bb[0],&SQ_ERI_F_D_F_D_b[0],&SQ_ERI_F_S_F_D[0],&SQ_ERI_F_F_G_D_bc[0],&SQ_ERI_F_F_D_D_b[0],&SQ_ERI_F_P_G_D_c[0],&SQ_ERI_F_P_D_D[0],abcd);

  hgp_os_eri_f_d_f_d_d2_deriv_4(&SQ_ERI_F_F_G_D_bc[0],&SQ_ERI_F_F_D_D_b[0],&SQ_ERI_F_P_G_D_c[0],&SQ_ERI_F_P_D_D[0],&SQ_ERI_F_D_H_D_cc[0],&SQ_ERI_F_D_F_D_c[0],&SQ_ERI_F_D_P_D[0],abcd);


}
