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
// BRA2
// X
// Y
// Z
// ####

void hgp_os_nai_f_sp_d1(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_NAI_G4x_S_C3_a = 0.0E0;
  Double I_NAI_G3xy_S_C3_a = 0.0E0;
  Double I_NAI_G3xz_S_C3_a = 0.0E0;
  Double I_NAI_G2x2y_S_C3_a = 0.0E0;
  Double I_NAI_G2xyz_S_C3_a = 0.0E0;
  Double I_NAI_G2x2z_S_C3_a = 0.0E0;
  Double I_NAI_Gx3y_S_C3_a = 0.0E0;
  Double I_NAI_Gx2yz_S_C3_a = 0.0E0;
  Double I_NAI_Gxy2z_S_C3_a = 0.0E0;
  Double I_NAI_Gx3z_S_C3_a = 0.0E0;
  Double I_NAI_G4y_S_C3_a = 0.0E0;
  Double I_NAI_G3yz_S_C3_a = 0.0E0;
  Double I_NAI_G2y2z_S_C3_a = 0.0E0;
  Double I_NAI_Gy3z_S_C3_a = 0.0E0;
  Double I_NAI_G4z_S_C3_a = 0.0E0;
  Double I_NAI_D2x_S_C3 = 0.0E0;
  Double I_NAI_Dxy_S_C3 = 0.0E0;
  Double I_NAI_Dxz_S_C3 = 0.0E0;
  Double I_NAI_D2y_S_C3 = 0.0E0;
  Double I_NAI_Dyz_S_C3 = 0.0E0;
  Double I_NAI_D2z_S_C3 = 0.0E0;
  Double I_NAI_F3x_S_C1003 = 0.0E0;
  Double I_NAI_F2xy_S_C1003 = 0.0E0;
  Double I_NAI_F2xz_S_C1003 = 0.0E0;
  Double I_NAI_Fx2y_S_C1003 = 0.0E0;
  Double I_NAI_Fxyz_S_C1003 = 0.0E0;
  Double I_NAI_Fx2z_S_C1003 = 0.0E0;
  Double I_NAI_F3y_S_C1003 = 0.0E0;
  Double I_NAI_F2yz_S_C1003 = 0.0E0;
  Double I_NAI_Fy2z_S_C1003 = 0.0E0;
  Double I_NAI_F3z_S_C1003 = 0.0E0;
  Double I_NAI_H5x_S_C1003_a = 0.0E0;
  Double I_NAI_H4xy_S_C1003_a = 0.0E0;
  Double I_NAI_H4xz_S_C1003_a = 0.0E0;
  Double I_NAI_H3x2y_S_C1003_a = 0.0E0;
  Double I_NAI_H3xyz_S_C1003_a = 0.0E0;
  Double I_NAI_H3x2z_S_C1003_a = 0.0E0;
  Double I_NAI_H2x3y_S_C1003_a = 0.0E0;
  Double I_NAI_H2x2yz_S_C1003_a = 0.0E0;
  Double I_NAI_H2xy2z_S_C1003_a = 0.0E0;
  Double I_NAI_H2x3z_S_C1003_a = 0.0E0;
  Double I_NAI_Hx4y_S_C1003_a = 0.0E0;
  Double I_NAI_Hx3yz_S_C1003_a = 0.0E0;
  Double I_NAI_Hx2y2z_S_C1003_a = 0.0E0;
  Double I_NAI_Hxy3z_S_C1003_a = 0.0E0;
  Double I_NAI_Hx4z_S_C1003_a = 0.0E0;
  Double I_NAI_H5y_S_C1003_a = 0.0E0;
  Double I_NAI_H4yz_S_C1003_a = 0.0E0;
  Double I_NAI_H3y2z_S_C1003_a = 0.0E0;
  Double I_NAI_H2y3z_S_C1003_a = 0.0E0;
  Double I_NAI_Hy4z_S_C1003_a = 0.0E0;
  Double I_NAI_H5z_S_C1003_a = 0.0E0;
  Double I_NAI_G4x_S_C1003_a = 0.0E0;
  Double I_NAI_G3xy_S_C1003_a = 0.0E0;
  Double I_NAI_G3xz_S_C1003_a = 0.0E0;
  Double I_NAI_G2x2y_S_C1003_a = 0.0E0;
  Double I_NAI_G2xyz_S_C1003_a = 0.0E0;
  Double I_NAI_G2x2z_S_C1003_a = 0.0E0;
  Double I_NAI_Gx3y_S_C1003_a = 0.0E0;
  Double I_NAI_Gx2yz_S_C1003_a = 0.0E0;
  Double I_NAI_Gxy2z_S_C1003_a = 0.0E0;
  Double I_NAI_Gx3z_S_C1003_a = 0.0E0;
  Double I_NAI_G4y_S_C1003_a = 0.0E0;
  Double I_NAI_G3yz_S_C1003_a = 0.0E0;
  Double I_NAI_G2y2z_S_C1003_a = 0.0E0;
  Double I_NAI_Gy3z_S_C1003_a = 0.0E0;
  Double I_NAI_G4z_S_C1003_a = 0.0E0;
  Double I_NAI_D2x_S_C1003 = 0.0E0;
  Double I_NAI_Dxy_S_C1003 = 0.0E0;
  Double I_NAI_Dxz_S_C1003 = 0.0E0;
  Double I_NAI_D2y_S_C1003 = 0.0E0;
  Double I_NAI_Dyz_S_C1003 = 0.0E0;
  Double I_NAI_D2z_S_C1003 = 0.0E0;
  Double I_NAI_G4x_S_C3_b = 0.0E0;
  Double I_NAI_G3xy_S_C3_b = 0.0E0;
  Double I_NAI_G3xz_S_C3_b = 0.0E0;
  Double I_NAI_G2x2y_S_C3_b = 0.0E0;
  Double I_NAI_G2xyz_S_C3_b = 0.0E0;
  Double I_NAI_G2x2z_S_C3_b = 0.0E0;
  Double I_NAI_Gx3y_S_C3_b = 0.0E0;
  Double I_NAI_Gx2yz_S_C3_b = 0.0E0;
  Double I_NAI_Gxy2z_S_C3_b = 0.0E0;
  Double I_NAI_Gx3z_S_C3_b = 0.0E0;
  Double I_NAI_G4y_S_C3_b = 0.0E0;
  Double I_NAI_G3yz_S_C3_b = 0.0E0;
  Double I_NAI_G2y2z_S_C3_b = 0.0E0;
  Double I_NAI_Gy3z_S_C3_b = 0.0E0;
  Double I_NAI_G4z_S_C3_b = 0.0E0;
  Double I_NAI_F3x_S_C3_b = 0.0E0;
  Double I_NAI_F2xy_S_C3_b = 0.0E0;
  Double I_NAI_F2xz_S_C3_b = 0.0E0;
  Double I_NAI_Fx2y_S_C3_b = 0.0E0;
  Double I_NAI_Fxyz_S_C3_b = 0.0E0;
  Double I_NAI_Fx2z_S_C3_b = 0.0E0;
  Double I_NAI_F3y_S_C3_b = 0.0E0;
  Double I_NAI_F2yz_S_C3_b = 0.0E0;
  Double I_NAI_Fy2z_S_C3_b = 0.0E0;
  Double I_NAI_F3z_S_C3_b = 0.0E0;
  Double I_NAI_H5x_S_C1003_b = 0.0E0;
  Double I_NAI_H4xy_S_C1003_b = 0.0E0;
  Double I_NAI_H4xz_S_C1003_b = 0.0E0;
  Double I_NAI_H3x2y_S_C1003_b = 0.0E0;
  Double I_NAI_H3xyz_S_C1003_b = 0.0E0;
  Double I_NAI_H3x2z_S_C1003_b = 0.0E0;
  Double I_NAI_H2x3y_S_C1003_b = 0.0E0;
  Double I_NAI_H2x2yz_S_C1003_b = 0.0E0;
  Double I_NAI_H2xy2z_S_C1003_b = 0.0E0;
  Double I_NAI_H2x3z_S_C1003_b = 0.0E0;
  Double I_NAI_Hx4y_S_C1003_b = 0.0E0;
  Double I_NAI_Hx3yz_S_C1003_b = 0.0E0;
  Double I_NAI_Hx2y2z_S_C1003_b = 0.0E0;
  Double I_NAI_Hxy3z_S_C1003_b = 0.0E0;
  Double I_NAI_Hx4z_S_C1003_b = 0.0E0;
  Double I_NAI_H5y_S_C1003_b = 0.0E0;
  Double I_NAI_H4yz_S_C1003_b = 0.0E0;
  Double I_NAI_H3y2z_S_C1003_b = 0.0E0;
  Double I_NAI_H2y3z_S_C1003_b = 0.0E0;
  Double I_NAI_Hy4z_S_C1003_b = 0.0E0;
  Double I_NAI_H5z_S_C1003_b = 0.0E0;
  Double I_NAI_G4x_S_C1003_b = 0.0E0;
  Double I_NAI_G3xy_S_C1003_b = 0.0E0;
  Double I_NAI_G3xz_S_C1003_b = 0.0E0;
  Double I_NAI_G2x2y_S_C1003_b = 0.0E0;
  Double I_NAI_G2xyz_S_C1003_b = 0.0E0;
  Double I_NAI_G2x2z_S_C1003_b = 0.0E0;
  Double I_NAI_Gx3y_S_C1003_b = 0.0E0;
  Double I_NAI_Gx2yz_S_C1003_b = 0.0E0;
  Double I_NAI_Gxy2z_S_C1003_b = 0.0E0;
  Double I_NAI_Gx3z_S_C1003_b = 0.0E0;
  Double I_NAI_G4y_S_C1003_b = 0.0E0;
  Double I_NAI_G3yz_S_C1003_b = 0.0E0;
  Double I_NAI_G2y2z_S_C1003_b = 0.0E0;
  Double I_NAI_Gy3z_S_C1003_b = 0.0E0;
  Double I_NAI_G4z_S_C1003_b = 0.0E0;
  Double I_NAI_F3x_S_C1003_b = 0.0E0;
  Double I_NAI_F2xy_S_C1003_b = 0.0E0;
  Double I_NAI_F2xz_S_C1003_b = 0.0E0;
  Double I_NAI_Fx2y_S_C1003_b = 0.0E0;
  Double I_NAI_Fxyz_S_C1003_b = 0.0E0;
  Double I_NAI_Fx2z_S_C1003_b = 0.0E0;
  Double I_NAI_F3y_S_C1003_b = 0.0E0;
  Double I_NAI_F2yz_S_C1003_b = 0.0E0;
  Double I_NAI_Fy2z_S_C1003_b = 0.0E0;
  Double I_NAI_F3z_S_C1003_b = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double ic2_1 = icoe[ip2+1*inp2];
    Double onedz = iexp[ip2];
    Double rho   = 1.0E0/onedz;
    Double zeta  = rho;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
    Double sqrho = sqrt(rho);
    Double fbra  = ifac[ip2];
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
      Double PNX   = PX - N[iAtom*3  ];
      Double PNY   = PY - N[iAtom*3+1];
      Double PNZ   = PZ - N[iAtom*3+2];
      Double PN2   = PNX*PNX+PNY*PNY+PNZ*PNZ;
      Double charge= Z[iAtom];
      Double u     = rho*PN2;
      Double squ   = sqrt(u);
      Double prefactor = -charge*fbra;

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

      Double I_NAI_S_S_vrr  = 0.0E0;
      Double I_NAI_S_S_M1_vrr  = 0.0E0;
      Double I_NAI_S_S_M2_vrr  = 0.0E0;
      Double I_NAI_S_S_M3_vrr  = 0.0E0;
      Double I_NAI_S_S_M4_vrr  = 0.0E0;
      Double I_NAI_S_S_M5_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_NAI_S_S_vrr_d  = 0.0E0;
      double I_NAI_S_S_M1_vrr_d  = 0.0E0;
      double I_NAI_S_S_M2_vrr_d  = 0.0E0;
      double I_NAI_S_S_M3_vrr_d  = 0.0E0;
      double I_NAI_S_S_M4_vrr_d  = 0.0E0;
      double I_NAI_S_S_M5_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER45;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER43*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER41*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER39*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER37*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER35*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER33*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER31*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER29*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER27*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER25*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER23*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER21*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER19*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER17*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER15*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = 1.0E0+u2*ONEOVER13*I_NAI_S_S_M5_vrr;
        I_NAI_S_S_M5_vrr = ONEOVER11*I_NAI_S_S_M5_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M5_vrr  = f*I_NAI_S_S_M5_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_NAI_S_S_M4_vrr  = ONEOVER9*(u2*I_NAI_S_S_M5_vrr+f);
        I_NAI_S_S_M3_vrr  = ONEOVER7*(u2*I_NAI_S_S_M4_vrr+f);
        I_NAI_S_S_M2_vrr  = ONEOVER5*(u2*I_NAI_S_S_M3_vrr+f);
        I_NAI_S_S_M1_vrr  = ONEOVER3*(u2*I_NAI_S_S_M2_vrr+f);
        I_NAI_S_S_vrr  = ONEOVER1*(u2*I_NAI_S_S_M1_vrr+f);

      }else{
#ifdef WITH_SINGLE_PRECISION

        // recompute the variable in terms of double accuracy
        double u_d     = u;
        double rho_d   = rho;
        double fac_d   = prefactor;
        double sqrho_d = sqrt(rho_d);
        double squ_d   = sqrt(u_d);

        // use erf function to get (SS|SS)^{0}
        if (fabs(u_d)<THRESHOLD_MATH) {
          I_NAI_S_S_vrr_d = fac_d*sqrho_d*TWOOVERSQRTPI;
        }else{
          I_NAI_S_S_vrr_d = (fac_d*sqrho_d/squ_d)*erf(squ_d);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        double oneO2u  = 0.5E0/u_d;
        double eu      = exp(-u_d);
        double f       = TWOOVERSQRTPI*fac_d*sqrho_d*eu;
        I_NAI_S_S_M1_vrr_d = oneO2u*(1.0E0*I_NAI_S_S_vrr_d-f);
        I_NAI_S_S_M2_vrr_d = oneO2u*(3.0E0*I_NAI_S_S_M1_vrr_d-f);
        I_NAI_S_S_M3_vrr_d = oneO2u*(5.0E0*I_NAI_S_S_M2_vrr_d-f);
        I_NAI_S_S_M4_vrr_d = oneO2u*(7.0E0*I_NAI_S_S_M3_vrr_d-f);
        I_NAI_S_S_M5_vrr_d = oneO2u*(9.0E0*I_NAI_S_S_M4_vrr_d-f);

        // write the double result back to the float var
        I_NAI_S_S_vrr = static_cast<Double>(I_NAI_S_S_vrr_d);
        I_NAI_S_S_M1_vrr = static_cast<Double>(I_NAI_S_S_M1_vrr_d);
        I_NAI_S_S_M2_vrr = static_cast<Double>(I_NAI_S_S_M2_vrr_d);
        I_NAI_S_S_M3_vrr = static_cast<Double>(I_NAI_S_S_M3_vrr_d);
        I_NAI_S_S_M4_vrr = static_cast<Double>(I_NAI_S_S_M4_vrr_d);
        I_NAI_S_S_M5_vrr = static_cast<Double>(I_NAI_S_S_M5_vrr_d);

#else

        // use erf function to get (SS|SS)^{0}
        if (fabs(u)<THRESHOLD_MATH) {
          I_NAI_S_S_vrr = prefactor*sqrho*TWOOVERSQRTPI;
        }else{
          I_NAI_S_S_vrr = (prefactor*sqrho/squ)*erf(squ);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        Double oneO2u = 0.5E0/u;
        Double eu     = exp(-u);
        Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M1_vrr = oneO2u*(1.0E0*I_NAI_S_S_vrr-f);
        I_NAI_S_S_M2_vrr = oneO2u*(3.0E0*I_NAI_S_S_M1_vrr-f);
        I_NAI_S_S_M3_vrr = oneO2u*(5.0E0*I_NAI_S_S_M2_vrr-f);
        I_NAI_S_S_M4_vrr = oneO2u*(7.0E0*I_NAI_S_S_M3_vrr-f);
        I_NAI_S_S_M5_vrr = oneO2u*(9.0E0*I_NAI_S_S_M4_vrr-f);

#endif

      }


      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M4
       * RHS shell quartet name: SQ_NAI_S_S_M5
       ************************************************************/
      Double I_NAI_Px_S_M4_vrr = PAX*I_NAI_S_S_M4_vrr-PNX*I_NAI_S_S_M5_vrr;
      Double I_NAI_Py_S_M4_vrr = PAY*I_NAI_S_S_M4_vrr-PNY*I_NAI_S_S_M5_vrr;
      Double I_NAI_Pz_S_M4_vrr = PAZ*I_NAI_S_S_M4_vrr-PNZ*I_NAI_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M3
       * RHS shell quartet name: SQ_NAI_S_S_M4
       ************************************************************/
      Double I_NAI_Px_S_M3_vrr = PAX*I_NAI_S_S_M3_vrr-PNX*I_NAI_S_S_M4_vrr;
      Double I_NAI_Py_S_M3_vrr = PAY*I_NAI_S_S_M3_vrr-PNY*I_NAI_S_S_M4_vrr;
      Double I_NAI_Pz_S_M3_vrr = PAZ*I_NAI_S_S_M3_vrr-PNZ*I_NAI_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M3
       * RHS shell quartet name: SQ_NAI_P_S_M4
       * RHS shell quartet name: SQ_NAI_S_S_M3
       * RHS shell quartet name: SQ_NAI_S_S_M4
       ************************************************************/
      Double I_NAI_D2x_S_M3_vrr = PAX*I_NAI_Px_S_M3_vrr-PNX*I_NAI_Px_S_M4_vrr+oned2z*I_NAI_S_S_M3_vrr-oned2z*I_NAI_S_S_M4_vrr;
      Double I_NAI_D2y_S_M3_vrr = PAY*I_NAI_Py_S_M3_vrr-PNY*I_NAI_Py_S_M4_vrr+oned2z*I_NAI_S_S_M3_vrr-oned2z*I_NAI_S_S_M4_vrr;
      Double I_NAI_D2z_S_M3_vrr = PAZ*I_NAI_Pz_S_M3_vrr-PNZ*I_NAI_Pz_S_M4_vrr+oned2z*I_NAI_S_S_M3_vrr-oned2z*I_NAI_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M3
       ************************************************************/
      Double I_NAI_Px_S_M2_vrr = PAX*I_NAI_S_S_M2_vrr-PNX*I_NAI_S_S_M3_vrr;
      Double I_NAI_Py_S_M2_vrr = PAY*I_NAI_S_S_M2_vrr-PNY*I_NAI_S_S_M3_vrr;
      Double I_NAI_Pz_S_M2_vrr = PAZ*I_NAI_S_S_M2_vrr-PNZ*I_NAI_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_P_S_M3
       * RHS shell quartet name: SQ_NAI_S_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M3
       ************************************************************/
      Double I_NAI_D2x_S_M2_vrr = PAX*I_NAI_Px_S_M2_vrr-PNX*I_NAI_Px_S_M3_vrr+oned2z*I_NAI_S_S_M2_vrr-oned2z*I_NAI_S_S_M3_vrr;
      Double I_NAI_D2y_S_M2_vrr = PAY*I_NAI_Py_S_M2_vrr-PNY*I_NAI_Py_S_M3_vrr+oned2z*I_NAI_S_S_M2_vrr-oned2z*I_NAI_S_S_M3_vrr;
      Double I_NAI_D2z_S_M2_vrr = PAZ*I_NAI_Pz_S_M2_vrr-PNZ*I_NAI_Pz_S_M3_vrr+oned2z*I_NAI_S_S_M2_vrr-oned2z*I_NAI_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M2
       * RHS shell quartet name: SQ_NAI_D_S_M3
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_P_S_M3
       ************************************************************/
      Double I_NAI_F3x_S_M2_vrr = PAX*I_NAI_D2x_S_M2_vrr-PNX*I_NAI_D2x_S_M3_vrr+2*oned2z*I_NAI_Px_S_M2_vrr-2*oned2z*I_NAI_Px_S_M3_vrr;
      Double I_NAI_F2xy_S_M2_vrr = PAY*I_NAI_D2x_S_M2_vrr-PNY*I_NAI_D2x_S_M3_vrr;
      Double I_NAI_F2xz_S_M2_vrr = PAZ*I_NAI_D2x_S_M2_vrr-PNZ*I_NAI_D2x_S_M3_vrr;
      Double I_NAI_F3y_S_M2_vrr = PAY*I_NAI_D2y_S_M2_vrr-PNY*I_NAI_D2y_S_M3_vrr+2*oned2z*I_NAI_Py_S_M2_vrr-2*oned2z*I_NAI_Py_S_M3_vrr;
      Double I_NAI_F2yz_S_M2_vrr = PAZ*I_NAI_D2y_S_M2_vrr-PNZ*I_NAI_D2y_S_M3_vrr;
      Double I_NAI_F3z_S_M2_vrr = PAZ*I_NAI_D2z_S_M2_vrr-PNZ*I_NAI_D2z_S_M3_vrr+2*oned2z*I_NAI_Pz_S_M2_vrr-2*oned2z*I_NAI_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M1
       * RHS shell quartet name: SQ_NAI_S_S_M2
       ************************************************************/
      Double I_NAI_Px_S_M1_vrr = PAX*I_NAI_S_S_M1_vrr-PNX*I_NAI_S_S_M2_vrr;
      Double I_NAI_Py_S_M1_vrr = PAY*I_NAI_S_S_M1_vrr-PNY*I_NAI_S_S_M2_vrr;
      Double I_NAI_Pz_S_M1_vrr = PAZ*I_NAI_S_S_M1_vrr-PNZ*I_NAI_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M1
       * RHS shell quartet name: SQ_NAI_S_S_M2
       ************************************************************/
      Double I_NAI_D2x_S_M1_vrr = PAX*I_NAI_Px_S_M1_vrr-PNX*I_NAI_Px_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
      Double I_NAI_Dxy_S_M1_vrr = PAY*I_NAI_Px_S_M1_vrr-PNY*I_NAI_Px_S_M2_vrr;
      Double I_NAI_D2y_S_M1_vrr = PAY*I_NAI_Py_S_M1_vrr-PNY*I_NAI_Py_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
      Double I_NAI_D2z_S_M1_vrr = PAZ*I_NAI_Pz_S_M1_vrr-PNZ*I_NAI_Pz_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_D_S_M2
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_P_S_M2
       ************************************************************/
      Double I_NAI_F3x_S_M1_vrr = PAX*I_NAI_D2x_S_M1_vrr-PNX*I_NAI_D2x_S_M2_vrr+2*oned2z*I_NAI_Px_S_M1_vrr-2*oned2z*I_NAI_Px_S_M2_vrr;
      Double I_NAI_F2xy_S_M1_vrr = PAY*I_NAI_D2x_S_M1_vrr-PNY*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_F2xz_S_M1_vrr = PAZ*I_NAI_D2x_S_M1_vrr-PNZ*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_Fx2y_S_M1_vrr = PAX*I_NAI_D2y_S_M1_vrr-PNX*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_Fx2z_S_M1_vrr = PAX*I_NAI_D2z_S_M1_vrr-PNX*I_NAI_D2z_S_M2_vrr;
      Double I_NAI_F3y_S_M1_vrr = PAY*I_NAI_D2y_S_M1_vrr-PNY*I_NAI_D2y_S_M2_vrr+2*oned2z*I_NAI_Py_S_M1_vrr-2*oned2z*I_NAI_Py_S_M2_vrr;
      Double I_NAI_F2yz_S_M1_vrr = PAZ*I_NAI_D2y_S_M1_vrr-PNZ*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_F3z_S_M1_vrr = PAZ*I_NAI_D2z_S_M1_vrr-PNZ*I_NAI_D2z_S_M2_vrr+2*oned2z*I_NAI_Pz_S_M1_vrr-2*oned2z*I_NAI_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S_M1
       * RHS shell quartet name: SQ_NAI_F_S_M2
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_D_S_M2
       ************************************************************/
      Double I_NAI_G4x_S_M1_vrr = PAX*I_NAI_F3x_S_M1_vrr-PNX*I_NAI_F3x_S_M2_vrr+3*oned2z*I_NAI_D2x_S_M1_vrr-3*oned2z*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_G3xy_S_M1_vrr = PAY*I_NAI_F3x_S_M1_vrr-PNY*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_G3xz_S_M1_vrr = PAZ*I_NAI_F3x_S_M1_vrr-PNZ*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_G2x2y_S_M1_vrr = PAY*I_NAI_F2xy_S_M1_vrr-PNY*I_NAI_F2xy_S_M2_vrr+oned2z*I_NAI_D2x_S_M1_vrr-oned2z*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_G2x2z_S_M1_vrr = PAZ*I_NAI_F2xz_S_M1_vrr-PNZ*I_NAI_F2xz_S_M2_vrr+oned2z*I_NAI_D2x_S_M1_vrr-oned2z*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_Gx3y_S_M1_vrr = PAX*I_NAI_F3y_S_M1_vrr-PNX*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_Gx3z_S_M1_vrr = PAX*I_NAI_F3z_S_M1_vrr-PNX*I_NAI_F3z_S_M2_vrr;
      Double I_NAI_G4y_S_M1_vrr = PAY*I_NAI_F3y_S_M1_vrr-PNY*I_NAI_F3y_S_M2_vrr+3*oned2z*I_NAI_D2y_S_M1_vrr-3*oned2z*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_G3yz_S_M1_vrr = PAZ*I_NAI_F3y_S_M1_vrr-PNZ*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_G2y2z_S_M1_vrr = PAZ*I_NAI_F2yz_S_M1_vrr-PNZ*I_NAI_F2yz_S_M2_vrr+oned2z*I_NAI_D2y_S_M1_vrr-oned2z*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_Gy3z_S_M1_vrr = PAY*I_NAI_F3z_S_M1_vrr-PNY*I_NAI_F3z_S_M2_vrr;
      Double I_NAI_G4z_S_M1_vrr = PAZ*I_NAI_F3z_S_M1_vrr-PNZ*I_NAI_F3z_S_M2_vrr+3*oned2z*I_NAI_D2z_S_M1_vrr-3*oned2z*I_NAI_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_Px_S_vrr = PAX*I_NAI_S_S_vrr-PNX*I_NAI_S_S_M1_vrr;
      Double I_NAI_Py_S_vrr = PAY*I_NAI_S_S_vrr-PNY*I_NAI_S_S_M1_vrr;
      Double I_NAI_Pz_S_vrr = PAZ*I_NAI_S_S_vrr-PNZ*I_NAI_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_D2x_S_vrr = PAX*I_NAI_Px_S_vrr-PNX*I_NAI_Px_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_Dxy_S_vrr = PAY*I_NAI_Px_S_vrr-PNY*I_NAI_Px_S_M1_vrr;
      Double I_NAI_Dxz_S_vrr = PAZ*I_NAI_Px_S_vrr-PNZ*I_NAI_Px_S_M1_vrr;
      Double I_NAI_D2y_S_vrr = PAY*I_NAI_Py_S_vrr-PNY*I_NAI_Py_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_Dyz_S_vrr = PAZ*I_NAI_Py_S_vrr-PNZ*I_NAI_Py_S_M1_vrr;
      Double I_NAI_D2z_S_vrr = PAZ*I_NAI_Pz_S_vrr-PNZ*I_NAI_Pz_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       ************************************************************/
      Double I_NAI_F3x_S_vrr = PAX*I_NAI_D2x_S_vrr-PNX*I_NAI_D2x_S_M1_vrr+2*oned2z*I_NAI_Px_S_vrr-2*oned2z*I_NAI_Px_S_M1_vrr;
      Double I_NAI_F2xy_S_vrr = PAY*I_NAI_D2x_S_vrr-PNY*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_F2xz_S_vrr = PAZ*I_NAI_D2x_S_vrr-PNZ*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_Fx2y_S_vrr = PAX*I_NAI_D2y_S_vrr-PNX*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_Fxyz_S_vrr = PAZ*I_NAI_Dxy_S_vrr-PNZ*I_NAI_Dxy_S_M1_vrr;
      Double I_NAI_Fx2z_S_vrr = PAX*I_NAI_D2z_S_vrr-PNX*I_NAI_D2z_S_M1_vrr;
      Double I_NAI_F3y_S_vrr = PAY*I_NAI_D2y_S_vrr-PNY*I_NAI_D2y_S_M1_vrr+2*oned2z*I_NAI_Py_S_vrr-2*oned2z*I_NAI_Py_S_M1_vrr;
      Double I_NAI_F2yz_S_vrr = PAZ*I_NAI_D2y_S_vrr-PNZ*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_Fy2z_S_vrr = PAY*I_NAI_D2z_S_vrr-PNY*I_NAI_D2z_S_M1_vrr;
      Double I_NAI_F3z_S_vrr = PAZ*I_NAI_D2z_S_vrr-PNZ*I_NAI_D2z_S_M1_vrr+2*oned2z*I_NAI_Pz_S_vrr-2*oned2z*I_NAI_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S
       * RHS shell quartet name: SQ_NAI_F_S_M1
       * RHS shell quartet name: SQ_NAI_D_S
       * RHS shell quartet name: SQ_NAI_D_S_M1
       ************************************************************/
      Double I_NAI_G4x_S_vrr = PAX*I_NAI_F3x_S_vrr-PNX*I_NAI_F3x_S_M1_vrr+3*oned2z*I_NAI_D2x_S_vrr-3*oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_G3xy_S_vrr = PAY*I_NAI_F3x_S_vrr-PNY*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_G3xz_S_vrr = PAZ*I_NAI_F3x_S_vrr-PNZ*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_G2x2y_S_vrr = PAY*I_NAI_F2xy_S_vrr-PNY*I_NAI_F2xy_S_M1_vrr+oned2z*I_NAI_D2x_S_vrr-oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_G2xyz_S_vrr = PAZ*I_NAI_F2xy_S_vrr-PNZ*I_NAI_F2xy_S_M1_vrr;
      Double I_NAI_G2x2z_S_vrr = PAZ*I_NAI_F2xz_S_vrr-PNZ*I_NAI_F2xz_S_M1_vrr+oned2z*I_NAI_D2x_S_vrr-oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_Gx3y_S_vrr = PAX*I_NAI_F3y_S_vrr-PNX*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_Gx2yz_S_vrr = PAZ*I_NAI_Fx2y_S_vrr-PNZ*I_NAI_Fx2y_S_M1_vrr;
      Double I_NAI_Gxy2z_S_vrr = PAY*I_NAI_Fx2z_S_vrr-PNY*I_NAI_Fx2z_S_M1_vrr;
      Double I_NAI_Gx3z_S_vrr = PAX*I_NAI_F3z_S_vrr-PNX*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_G4y_S_vrr = PAY*I_NAI_F3y_S_vrr-PNY*I_NAI_F3y_S_M1_vrr+3*oned2z*I_NAI_D2y_S_vrr-3*oned2z*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_G3yz_S_vrr = PAZ*I_NAI_F3y_S_vrr-PNZ*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_G2y2z_S_vrr = PAZ*I_NAI_F2yz_S_vrr-PNZ*I_NAI_F2yz_S_M1_vrr+oned2z*I_NAI_D2y_S_vrr-oned2z*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_Gy3z_S_vrr = PAY*I_NAI_F3z_S_vrr-PNY*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_G4z_S_vrr = PAZ*I_NAI_F3z_S_vrr-PNZ*I_NAI_F3z_S_M1_vrr+3*oned2z*I_NAI_D2z_S_vrr-3*oned2z*I_NAI_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_G_S
       * RHS shell quartet name: SQ_NAI_G_S_M1
       * RHS shell quartet name: SQ_NAI_F_S
       * RHS shell quartet name: SQ_NAI_F_S_M1
       ************************************************************/
      Double I_NAI_H5x_S_vrr = PAX*I_NAI_G4x_S_vrr-PNX*I_NAI_G4x_S_M1_vrr+4*oned2z*I_NAI_F3x_S_vrr-4*oned2z*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_H4xy_S_vrr = PAY*I_NAI_G4x_S_vrr-PNY*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_H4xz_S_vrr = PAZ*I_NAI_G4x_S_vrr-PNZ*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_H3x2y_S_vrr = PAY*I_NAI_G3xy_S_vrr-PNY*I_NAI_G3xy_S_M1_vrr+oned2z*I_NAI_F3x_S_vrr-oned2z*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_H3xyz_S_vrr = PAZ*I_NAI_G3xy_S_vrr-PNZ*I_NAI_G3xy_S_M1_vrr;
      Double I_NAI_H3x2z_S_vrr = PAZ*I_NAI_G3xz_S_vrr-PNZ*I_NAI_G3xz_S_M1_vrr+oned2z*I_NAI_F3x_S_vrr-oned2z*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_H2x3y_S_vrr = PAX*I_NAI_Gx3y_S_vrr-PNX*I_NAI_Gx3y_S_M1_vrr+oned2z*I_NAI_F3y_S_vrr-oned2z*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_H2x2yz_S_vrr = PAZ*I_NAI_G2x2y_S_vrr-PNZ*I_NAI_G2x2y_S_M1_vrr;
      Double I_NAI_H2xy2z_S_vrr = PAY*I_NAI_G2x2z_S_vrr-PNY*I_NAI_G2x2z_S_M1_vrr;
      Double I_NAI_H2x3z_S_vrr = PAX*I_NAI_Gx3z_S_vrr-PNX*I_NAI_Gx3z_S_M1_vrr+oned2z*I_NAI_F3z_S_vrr-oned2z*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_Hx4y_S_vrr = PAX*I_NAI_G4y_S_vrr-PNX*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_Hx3yz_S_vrr = PAZ*I_NAI_Gx3y_S_vrr-PNZ*I_NAI_Gx3y_S_M1_vrr;
      Double I_NAI_Hx2y2z_S_vrr = PAX*I_NAI_G2y2z_S_vrr-PNX*I_NAI_G2y2z_S_M1_vrr;
      Double I_NAI_Hxy3z_S_vrr = PAY*I_NAI_Gx3z_S_vrr-PNY*I_NAI_Gx3z_S_M1_vrr;
      Double I_NAI_Hx4z_S_vrr = PAX*I_NAI_G4z_S_vrr-PNX*I_NAI_G4z_S_M1_vrr;
      Double I_NAI_H5y_S_vrr = PAY*I_NAI_G4y_S_vrr-PNY*I_NAI_G4y_S_M1_vrr+4*oned2z*I_NAI_F3y_S_vrr-4*oned2z*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_H4yz_S_vrr = PAZ*I_NAI_G4y_S_vrr-PNZ*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_H3y2z_S_vrr = PAZ*I_NAI_G3yz_S_vrr-PNZ*I_NAI_G3yz_S_M1_vrr+oned2z*I_NAI_F3y_S_vrr-oned2z*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_H2y3z_S_vrr = PAY*I_NAI_Gy3z_S_vrr-PNY*I_NAI_Gy3z_S_M1_vrr+oned2z*I_NAI_F3z_S_vrr-oned2z*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_Hy4z_S_vrr = PAY*I_NAI_G4z_S_vrr-PNY*I_NAI_G4z_S_M1_vrr;
      Double I_NAI_H5z_S_vrr = PAZ*I_NAI_G4z_S_vrr-PNZ*I_NAI_G4z_S_M1_vrr+4*oned2z*I_NAI_F3z_S_vrr-4*oned2z*I_NAI_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C3_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C3_a_coefs = ic2*alpha;
      I_NAI_G4x_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C3_a += SQ_NAI_G_S_C3_a_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C3
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C3_coefs = ic2;
      I_NAI_D2x_S_C3 += SQ_NAI_D_S_C3_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C3 += SQ_NAI_D_S_C3_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C3 += SQ_NAI_D_S_C3_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C3 += SQ_NAI_D_S_C3_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C3 += SQ_NAI_D_S_C3_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C3 += SQ_NAI_D_S_C3_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1003_coefs = ic2_1;
      I_NAI_F3x_S_C1003 += SQ_NAI_F_S_C1003_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1003 += SQ_NAI_F_S_C1003_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1003 += SQ_NAI_F_S_C1003_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1003 += SQ_NAI_F_S_C1003_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1003 += SQ_NAI_F_S_C1003_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1003 += SQ_NAI_F_S_C1003_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1003 += SQ_NAI_F_S_C1003_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1003 += SQ_NAI_F_S_C1003_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1003 += SQ_NAI_F_S_C1003_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1003 += SQ_NAI_F_S_C1003_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_C1003_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C1003_a_coefs = ic2_1*alpha;
      I_NAI_H5x_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C1003_a += SQ_NAI_H_S_C1003_a_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1003_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1003_a_coefs = ic2_1*alpha;
      I_NAI_G4x_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1003_a += SQ_NAI_G_S_C1003_a_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1003_coefs = ic2_1;
      I_NAI_D2x_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1003 += SQ_NAI_D_S_C1003_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C3_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C3_b_coefs = ic2*beta;
      I_NAI_G4x_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C3_b += SQ_NAI_G_S_C3_b_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C3_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C3_b_coefs = ic2*beta;
      I_NAI_F3x_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C3_b += SQ_NAI_F_S_C3_b_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_C1003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_C1003_b_coefs = ic2_1*beta;
      I_NAI_H5x_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_C1003_b += SQ_NAI_H_S_C1003_b_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1003_b_coefs = ic2_1*beta;
      I_NAI_G4x_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1003_b += SQ_NAI_G_S_C1003_b_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1003_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1003_b_coefs = ic2_1*beta;
      I_NAI_F3x_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1003_b += SQ_NAI_F_S_C1003_b_coefs*I_NAI_F3z_S_vrr;
    }
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
   * shell quartet name: SQ_NAI_D_P_C1003
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1003
   * RHS shell quartet name: SQ_NAI_D_S_C1003
   ************************************************************/
  Double I_NAI_D2x_Px_C1003 = I_NAI_F3x_S_C1003+ABX*I_NAI_D2x_S_C1003;
  Double I_NAI_Dxy_Px_C1003 = I_NAI_F2xy_S_C1003+ABX*I_NAI_Dxy_S_C1003;
  Double I_NAI_Dxz_Px_C1003 = I_NAI_F2xz_S_C1003+ABX*I_NAI_Dxz_S_C1003;
  Double I_NAI_D2y_Px_C1003 = I_NAI_Fx2y_S_C1003+ABX*I_NAI_D2y_S_C1003;
  Double I_NAI_Dyz_Px_C1003 = I_NAI_Fxyz_S_C1003+ABX*I_NAI_Dyz_S_C1003;
  Double I_NAI_D2z_Px_C1003 = I_NAI_Fx2z_S_C1003+ABX*I_NAI_D2z_S_C1003;
  Double I_NAI_D2x_Py_C1003 = I_NAI_F2xy_S_C1003+ABY*I_NAI_D2x_S_C1003;
  Double I_NAI_Dxy_Py_C1003 = I_NAI_Fx2y_S_C1003+ABY*I_NAI_Dxy_S_C1003;
  Double I_NAI_Dxz_Py_C1003 = I_NAI_Fxyz_S_C1003+ABY*I_NAI_Dxz_S_C1003;
  Double I_NAI_D2y_Py_C1003 = I_NAI_F3y_S_C1003+ABY*I_NAI_D2y_S_C1003;
  Double I_NAI_Dyz_Py_C1003 = I_NAI_F2yz_S_C1003+ABY*I_NAI_Dyz_S_C1003;
  Double I_NAI_D2z_Py_C1003 = I_NAI_Fy2z_S_C1003+ABY*I_NAI_D2z_S_C1003;
  Double I_NAI_D2x_Pz_C1003 = I_NAI_F2xz_S_C1003+ABZ*I_NAI_D2x_S_C1003;
  Double I_NAI_Dxy_Pz_C1003 = I_NAI_Fxyz_S_C1003+ABZ*I_NAI_Dxy_S_C1003;
  Double I_NAI_Dxz_Pz_C1003 = I_NAI_Fx2z_S_C1003+ABZ*I_NAI_Dxz_S_C1003;
  Double I_NAI_D2y_Pz_C1003 = I_NAI_F2yz_S_C1003+ABZ*I_NAI_D2y_S_C1003;
  Double I_NAI_Dyz_Pz_C1003 = I_NAI_Fy2z_S_C1003+ABZ*I_NAI_Dyz_S_C1003;
  Double I_NAI_D2z_Pz_C1003 = I_NAI_F3z_S_C1003+ABZ*I_NAI_D2z_S_C1003;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_C1003_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C1003_a
   * RHS shell quartet name: SQ_NAI_G_S_C1003_a
   ************************************************************/
  Double I_NAI_G4x_Px_C1003_a = I_NAI_H5x_S_C1003_a+ABX*I_NAI_G4x_S_C1003_a;
  Double I_NAI_G3xy_Px_C1003_a = I_NAI_H4xy_S_C1003_a+ABX*I_NAI_G3xy_S_C1003_a;
  Double I_NAI_G3xz_Px_C1003_a = I_NAI_H4xz_S_C1003_a+ABX*I_NAI_G3xz_S_C1003_a;
  Double I_NAI_G2x2y_Px_C1003_a = I_NAI_H3x2y_S_C1003_a+ABX*I_NAI_G2x2y_S_C1003_a;
  Double I_NAI_G2xyz_Px_C1003_a = I_NAI_H3xyz_S_C1003_a+ABX*I_NAI_G2xyz_S_C1003_a;
  Double I_NAI_G2x2z_Px_C1003_a = I_NAI_H3x2z_S_C1003_a+ABX*I_NAI_G2x2z_S_C1003_a;
  Double I_NAI_Gx3y_Px_C1003_a = I_NAI_H2x3y_S_C1003_a+ABX*I_NAI_Gx3y_S_C1003_a;
  Double I_NAI_Gx2yz_Px_C1003_a = I_NAI_H2x2yz_S_C1003_a+ABX*I_NAI_Gx2yz_S_C1003_a;
  Double I_NAI_Gxy2z_Px_C1003_a = I_NAI_H2xy2z_S_C1003_a+ABX*I_NAI_Gxy2z_S_C1003_a;
  Double I_NAI_Gx3z_Px_C1003_a = I_NAI_H2x3z_S_C1003_a+ABX*I_NAI_Gx3z_S_C1003_a;
  Double I_NAI_G4y_Px_C1003_a = I_NAI_Hx4y_S_C1003_a+ABX*I_NAI_G4y_S_C1003_a;
  Double I_NAI_G3yz_Px_C1003_a = I_NAI_Hx3yz_S_C1003_a+ABX*I_NAI_G3yz_S_C1003_a;
  Double I_NAI_G2y2z_Px_C1003_a = I_NAI_Hx2y2z_S_C1003_a+ABX*I_NAI_G2y2z_S_C1003_a;
  Double I_NAI_Gy3z_Px_C1003_a = I_NAI_Hxy3z_S_C1003_a+ABX*I_NAI_Gy3z_S_C1003_a;
  Double I_NAI_G4z_Px_C1003_a = I_NAI_Hx4z_S_C1003_a+ABX*I_NAI_G4z_S_C1003_a;
  Double I_NAI_G4x_Py_C1003_a = I_NAI_H4xy_S_C1003_a+ABY*I_NAI_G4x_S_C1003_a;
  Double I_NAI_G3xy_Py_C1003_a = I_NAI_H3x2y_S_C1003_a+ABY*I_NAI_G3xy_S_C1003_a;
  Double I_NAI_G3xz_Py_C1003_a = I_NAI_H3xyz_S_C1003_a+ABY*I_NAI_G3xz_S_C1003_a;
  Double I_NAI_G2x2y_Py_C1003_a = I_NAI_H2x3y_S_C1003_a+ABY*I_NAI_G2x2y_S_C1003_a;
  Double I_NAI_G2xyz_Py_C1003_a = I_NAI_H2x2yz_S_C1003_a+ABY*I_NAI_G2xyz_S_C1003_a;
  Double I_NAI_G2x2z_Py_C1003_a = I_NAI_H2xy2z_S_C1003_a+ABY*I_NAI_G2x2z_S_C1003_a;
  Double I_NAI_Gx3y_Py_C1003_a = I_NAI_Hx4y_S_C1003_a+ABY*I_NAI_Gx3y_S_C1003_a;
  Double I_NAI_Gx2yz_Py_C1003_a = I_NAI_Hx3yz_S_C1003_a+ABY*I_NAI_Gx2yz_S_C1003_a;
  Double I_NAI_Gxy2z_Py_C1003_a = I_NAI_Hx2y2z_S_C1003_a+ABY*I_NAI_Gxy2z_S_C1003_a;
  Double I_NAI_Gx3z_Py_C1003_a = I_NAI_Hxy3z_S_C1003_a+ABY*I_NAI_Gx3z_S_C1003_a;
  Double I_NAI_G4y_Py_C1003_a = I_NAI_H5y_S_C1003_a+ABY*I_NAI_G4y_S_C1003_a;
  Double I_NAI_G3yz_Py_C1003_a = I_NAI_H4yz_S_C1003_a+ABY*I_NAI_G3yz_S_C1003_a;
  Double I_NAI_G2y2z_Py_C1003_a = I_NAI_H3y2z_S_C1003_a+ABY*I_NAI_G2y2z_S_C1003_a;
  Double I_NAI_Gy3z_Py_C1003_a = I_NAI_H2y3z_S_C1003_a+ABY*I_NAI_Gy3z_S_C1003_a;
  Double I_NAI_G4z_Py_C1003_a = I_NAI_Hy4z_S_C1003_a+ABY*I_NAI_G4z_S_C1003_a;
  Double I_NAI_G4x_Pz_C1003_a = I_NAI_H4xz_S_C1003_a+ABZ*I_NAI_G4x_S_C1003_a;
  Double I_NAI_G3xy_Pz_C1003_a = I_NAI_H3xyz_S_C1003_a+ABZ*I_NAI_G3xy_S_C1003_a;
  Double I_NAI_G3xz_Pz_C1003_a = I_NAI_H3x2z_S_C1003_a+ABZ*I_NAI_G3xz_S_C1003_a;
  Double I_NAI_G2x2y_Pz_C1003_a = I_NAI_H2x2yz_S_C1003_a+ABZ*I_NAI_G2x2y_S_C1003_a;
  Double I_NAI_G2xyz_Pz_C1003_a = I_NAI_H2xy2z_S_C1003_a+ABZ*I_NAI_G2xyz_S_C1003_a;
  Double I_NAI_G2x2z_Pz_C1003_a = I_NAI_H2x3z_S_C1003_a+ABZ*I_NAI_G2x2z_S_C1003_a;
  Double I_NAI_Gx3y_Pz_C1003_a = I_NAI_Hx3yz_S_C1003_a+ABZ*I_NAI_Gx3y_S_C1003_a;
  Double I_NAI_Gx2yz_Pz_C1003_a = I_NAI_Hx2y2z_S_C1003_a+ABZ*I_NAI_Gx2yz_S_C1003_a;
  Double I_NAI_Gxy2z_Pz_C1003_a = I_NAI_Hxy3z_S_C1003_a+ABZ*I_NAI_Gxy2z_S_C1003_a;
  Double I_NAI_Gx3z_Pz_C1003_a = I_NAI_Hx4z_S_C1003_a+ABZ*I_NAI_Gx3z_S_C1003_a;
  Double I_NAI_G4y_Pz_C1003_a = I_NAI_H4yz_S_C1003_a+ABZ*I_NAI_G4y_S_C1003_a;
  Double I_NAI_G3yz_Pz_C1003_a = I_NAI_H3y2z_S_C1003_a+ABZ*I_NAI_G3yz_S_C1003_a;
  Double I_NAI_G2y2z_Pz_C1003_a = I_NAI_H2y3z_S_C1003_a+ABZ*I_NAI_G2y2z_S_C1003_a;
  Double I_NAI_Gy3z_Pz_C1003_a = I_NAI_Hy4z_S_C1003_a+ABZ*I_NAI_Gy3z_S_C1003_a;
  Double I_NAI_G4z_Pz_C1003_a = I_NAI_H5z_S_C1003_a+ABZ*I_NAI_G4z_S_C1003_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C3_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C3_b
   * RHS shell quartet name: SQ_NAI_F_S_C3_b
   ************************************************************/
  Double I_NAI_F3x_Px_C3_b = I_NAI_G4x_S_C3_b+ABX*I_NAI_F3x_S_C3_b;
  Double I_NAI_F2xy_Px_C3_b = I_NAI_G3xy_S_C3_b+ABX*I_NAI_F2xy_S_C3_b;
  Double I_NAI_F2xz_Px_C3_b = I_NAI_G3xz_S_C3_b+ABX*I_NAI_F2xz_S_C3_b;
  Double I_NAI_Fx2y_Px_C3_b = I_NAI_G2x2y_S_C3_b+ABX*I_NAI_Fx2y_S_C3_b;
  Double I_NAI_Fxyz_Px_C3_b = I_NAI_G2xyz_S_C3_b+ABX*I_NAI_Fxyz_S_C3_b;
  Double I_NAI_Fx2z_Px_C3_b = I_NAI_G2x2z_S_C3_b+ABX*I_NAI_Fx2z_S_C3_b;
  Double I_NAI_F3y_Px_C3_b = I_NAI_Gx3y_S_C3_b+ABX*I_NAI_F3y_S_C3_b;
  Double I_NAI_F2yz_Px_C3_b = I_NAI_Gx2yz_S_C3_b+ABX*I_NAI_F2yz_S_C3_b;
  Double I_NAI_Fy2z_Px_C3_b = I_NAI_Gxy2z_S_C3_b+ABX*I_NAI_Fy2z_S_C3_b;
  Double I_NAI_F3z_Px_C3_b = I_NAI_Gx3z_S_C3_b+ABX*I_NAI_F3z_S_C3_b;
  Double I_NAI_F3x_Py_C3_b = I_NAI_G3xy_S_C3_b+ABY*I_NAI_F3x_S_C3_b;
  Double I_NAI_F2xy_Py_C3_b = I_NAI_G2x2y_S_C3_b+ABY*I_NAI_F2xy_S_C3_b;
  Double I_NAI_F2xz_Py_C3_b = I_NAI_G2xyz_S_C3_b+ABY*I_NAI_F2xz_S_C3_b;
  Double I_NAI_Fx2y_Py_C3_b = I_NAI_Gx3y_S_C3_b+ABY*I_NAI_Fx2y_S_C3_b;
  Double I_NAI_Fxyz_Py_C3_b = I_NAI_Gx2yz_S_C3_b+ABY*I_NAI_Fxyz_S_C3_b;
  Double I_NAI_Fx2z_Py_C3_b = I_NAI_Gxy2z_S_C3_b+ABY*I_NAI_Fx2z_S_C3_b;
  Double I_NAI_F3y_Py_C3_b = I_NAI_G4y_S_C3_b+ABY*I_NAI_F3y_S_C3_b;
  Double I_NAI_F2yz_Py_C3_b = I_NAI_G3yz_S_C3_b+ABY*I_NAI_F2yz_S_C3_b;
  Double I_NAI_Fy2z_Py_C3_b = I_NAI_G2y2z_S_C3_b+ABY*I_NAI_Fy2z_S_C3_b;
  Double I_NAI_F3z_Py_C3_b = I_NAI_Gy3z_S_C3_b+ABY*I_NAI_F3z_S_C3_b;
  Double I_NAI_F3x_Pz_C3_b = I_NAI_G3xz_S_C3_b+ABZ*I_NAI_F3x_S_C3_b;
  Double I_NAI_F2xy_Pz_C3_b = I_NAI_G2xyz_S_C3_b+ABZ*I_NAI_F2xy_S_C3_b;
  Double I_NAI_F2xz_Pz_C3_b = I_NAI_G2x2z_S_C3_b+ABZ*I_NAI_F2xz_S_C3_b;
  Double I_NAI_Fx2y_Pz_C3_b = I_NAI_Gx2yz_S_C3_b+ABZ*I_NAI_Fx2y_S_C3_b;
  Double I_NAI_Fxyz_Pz_C3_b = I_NAI_Gxy2z_S_C3_b+ABZ*I_NAI_Fxyz_S_C3_b;
  Double I_NAI_Fx2z_Pz_C3_b = I_NAI_Gx3z_S_C3_b+ABZ*I_NAI_Fx2z_S_C3_b;
  Double I_NAI_F3y_Pz_C3_b = I_NAI_G3yz_S_C3_b+ABZ*I_NAI_F3y_S_C3_b;
  Double I_NAI_F2yz_Pz_C3_b = I_NAI_G2y2z_S_C3_b+ABZ*I_NAI_F2yz_S_C3_b;
  Double I_NAI_Fy2z_Pz_C3_b = I_NAI_Gy3z_S_C3_b+ABZ*I_NAI_Fy2z_S_C3_b;
  Double I_NAI_F3z_Pz_C3_b = I_NAI_G4z_S_C3_b+ABZ*I_NAI_F3z_S_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C1003_b
   * RHS shell quartet name: SQ_NAI_F_S_C1003_b
   ************************************************************/
  Double I_NAI_F3x_Px_C1003_b = I_NAI_G4x_S_C1003_b+ABX*I_NAI_F3x_S_C1003_b;
  Double I_NAI_F2xy_Px_C1003_b = I_NAI_G3xy_S_C1003_b+ABX*I_NAI_F2xy_S_C1003_b;
  Double I_NAI_F2xz_Px_C1003_b = I_NAI_G3xz_S_C1003_b+ABX*I_NAI_F2xz_S_C1003_b;
  Double I_NAI_Fx2y_Px_C1003_b = I_NAI_G2x2y_S_C1003_b+ABX*I_NAI_Fx2y_S_C1003_b;
  Double I_NAI_Fxyz_Px_C1003_b = I_NAI_G2xyz_S_C1003_b+ABX*I_NAI_Fxyz_S_C1003_b;
  Double I_NAI_Fx2z_Px_C1003_b = I_NAI_G2x2z_S_C1003_b+ABX*I_NAI_Fx2z_S_C1003_b;
  Double I_NAI_F3y_Px_C1003_b = I_NAI_Gx3y_S_C1003_b+ABX*I_NAI_F3y_S_C1003_b;
  Double I_NAI_F2yz_Px_C1003_b = I_NAI_Gx2yz_S_C1003_b+ABX*I_NAI_F2yz_S_C1003_b;
  Double I_NAI_Fy2z_Px_C1003_b = I_NAI_Gxy2z_S_C1003_b+ABX*I_NAI_Fy2z_S_C1003_b;
  Double I_NAI_F3z_Px_C1003_b = I_NAI_Gx3z_S_C1003_b+ABX*I_NAI_F3z_S_C1003_b;
  Double I_NAI_F3x_Py_C1003_b = I_NAI_G3xy_S_C1003_b+ABY*I_NAI_F3x_S_C1003_b;
  Double I_NAI_F2xy_Py_C1003_b = I_NAI_G2x2y_S_C1003_b+ABY*I_NAI_F2xy_S_C1003_b;
  Double I_NAI_F2xz_Py_C1003_b = I_NAI_G2xyz_S_C1003_b+ABY*I_NAI_F2xz_S_C1003_b;
  Double I_NAI_Fx2y_Py_C1003_b = I_NAI_Gx3y_S_C1003_b+ABY*I_NAI_Fx2y_S_C1003_b;
  Double I_NAI_Fxyz_Py_C1003_b = I_NAI_Gx2yz_S_C1003_b+ABY*I_NAI_Fxyz_S_C1003_b;
  Double I_NAI_Fx2z_Py_C1003_b = I_NAI_Gxy2z_S_C1003_b+ABY*I_NAI_Fx2z_S_C1003_b;
  Double I_NAI_F3y_Py_C1003_b = I_NAI_G4y_S_C1003_b+ABY*I_NAI_F3y_S_C1003_b;
  Double I_NAI_F2yz_Py_C1003_b = I_NAI_G3yz_S_C1003_b+ABY*I_NAI_F2yz_S_C1003_b;
  Double I_NAI_Fy2z_Py_C1003_b = I_NAI_G2y2z_S_C1003_b+ABY*I_NAI_Fy2z_S_C1003_b;
  Double I_NAI_F3z_Py_C1003_b = I_NAI_Gy3z_S_C1003_b+ABY*I_NAI_F3z_S_C1003_b;
  Double I_NAI_F3x_Pz_C1003_b = I_NAI_G3xz_S_C1003_b+ABZ*I_NAI_F3x_S_C1003_b;
  Double I_NAI_F2xy_Pz_C1003_b = I_NAI_G2xyz_S_C1003_b+ABZ*I_NAI_F2xy_S_C1003_b;
  Double I_NAI_F2xz_Pz_C1003_b = I_NAI_G2x2z_S_C1003_b+ABZ*I_NAI_F2xz_S_C1003_b;
  Double I_NAI_Fx2y_Pz_C1003_b = I_NAI_Gx2yz_S_C1003_b+ABZ*I_NAI_Fx2y_S_C1003_b;
  Double I_NAI_Fxyz_Pz_C1003_b = I_NAI_Gxy2z_S_C1003_b+ABZ*I_NAI_Fxyz_S_C1003_b;
  Double I_NAI_Fx2z_Pz_C1003_b = I_NAI_Gx3z_S_C1003_b+ABZ*I_NAI_Fx2z_S_C1003_b;
  Double I_NAI_F3y_Pz_C1003_b = I_NAI_G3yz_S_C1003_b+ABZ*I_NAI_F3y_S_C1003_b;
  Double I_NAI_F2yz_Pz_C1003_b = I_NAI_G2y2z_S_C1003_b+ABZ*I_NAI_F2yz_S_C1003_b;
  Double I_NAI_Fy2z_Pz_C1003_b = I_NAI_Gy3z_S_C1003_b+ABZ*I_NAI_Fy2z_S_C1003_b;
  Double I_NAI_F3z_Pz_C1003_b = I_NAI_G4z_S_C1003_b+ABZ*I_NAI_F3z_S_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_C1003_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 6 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_C1003_b
   * RHS shell quartet name: SQ_NAI_G_S_C1003_b
   ************************************************************/
  Double I_NAI_G4x_Px_C1003_b = I_NAI_H5x_S_C1003_b+ABX*I_NAI_G4x_S_C1003_b;
  Double I_NAI_G3xy_Px_C1003_b = I_NAI_H4xy_S_C1003_b+ABX*I_NAI_G3xy_S_C1003_b;
  Double I_NAI_G3xz_Px_C1003_b = I_NAI_H4xz_S_C1003_b+ABX*I_NAI_G3xz_S_C1003_b;
  Double I_NAI_G2x2y_Px_C1003_b = I_NAI_H3x2y_S_C1003_b+ABX*I_NAI_G2x2y_S_C1003_b;
  Double I_NAI_G2xyz_Px_C1003_b = I_NAI_H3xyz_S_C1003_b+ABX*I_NAI_G2xyz_S_C1003_b;
  Double I_NAI_G2x2z_Px_C1003_b = I_NAI_H3x2z_S_C1003_b+ABX*I_NAI_G2x2z_S_C1003_b;
  Double I_NAI_Gx3y_Px_C1003_b = I_NAI_H2x3y_S_C1003_b+ABX*I_NAI_Gx3y_S_C1003_b;
  Double I_NAI_Gx2yz_Px_C1003_b = I_NAI_H2x2yz_S_C1003_b+ABX*I_NAI_Gx2yz_S_C1003_b;
  Double I_NAI_Gxy2z_Px_C1003_b = I_NAI_H2xy2z_S_C1003_b+ABX*I_NAI_Gxy2z_S_C1003_b;
  Double I_NAI_Gx3z_Px_C1003_b = I_NAI_H2x3z_S_C1003_b+ABX*I_NAI_Gx3z_S_C1003_b;
  Double I_NAI_G4y_Px_C1003_b = I_NAI_Hx4y_S_C1003_b+ABX*I_NAI_G4y_S_C1003_b;
  Double I_NAI_G3yz_Px_C1003_b = I_NAI_Hx3yz_S_C1003_b+ABX*I_NAI_G3yz_S_C1003_b;
  Double I_NAI_G2y2z_Px_C1003_b = I_NAI_Hx2y2z_S_C1003_b+ABX*I_NAI_G2y2z_S_C1003_b;
  Double I_NAI_Gy3z_Px_C1003_b = I_NAI_Hxy3z_S_C1003_b+ABX*I_NAI_Gy3z_S_C1003_b;
  Double I_NAI_G4z_Px_C1003_b = I_NAI_Hx4z_S_C1003_b+ABX*I_NAI_G4z_S_C1003_b;
  Double I_NAI_G3xy_Py_C1003_b = I_NAI_H3x2y_S_C1003_b+ABY*I_NAI_G3xy_S_C1003_b;
  Double I_NAI_G3xz_Py_C1003_b = I_NAI_H3xyz_S_C1003_b+ABY*I_NAI_G3xz_S_C1003_b;
  Double I_NAI_G2x2y_Py_C1003_b = I_NAI_H2x3y_S_C1003_b+ABY*I_NAI_G2x2y_S_C1003_b;
  Double I_NAI_G2xyz_Py_C1003_b = I_NAI_H2x2yz_S_C1003_b+ABY*I_NAI_G2xyz_S_C1003_b;
  Double I_NAI_G2x2z_Py_C1003_b = I_NAI_H2xy2z_S_C1003_b+ABY*I_NAI_G2x2z_S_C1003_b;
  Double I_NAI_Gx3y_Py_C1003_b = I_NAI_Hx4y_S_C1003_b+ABY*I_NAI_Gx3y_S_C1003_b;
  Double I_NAI_Gx2yz_Py_C1003_b = I_NAI_Hx3yz_S_C1003_b+ABY*I_NAI_Gx2yz_S_C1003_b;
  Double I_NAI_Gxy2z_Py_C1003_b = I_NAI_Hx2y2z_S_C1003_b+ABY*I_NAI_Gxy2z_S_C1003_b;
  Double I_NAI_Gx3z_Py_C1003_b = I_NAI_Hxy3z_S_C1003_b+ABY*I_NAI_Gx3z_S_C1003_b;
  Double I_NAI_G4y_Py_C1003_b = I_NAI_H5y_S_C1003_b+ABY*I_NAI_G4y_S_C1003_b;
  Double I_NAI_G3yz_Py_C1003_b = I_NAI_H4yz_S_C1003_b+ABY*I_NAI_G3yz_S_C1003_b;
  Double I_NAI_G2y2z_Py_C1003_b = I_NAI_H3y2z_S_C1003_b+ABY*I_NAI_G2y2z_S_C1003_b;
  Double I_NAI_Gy3z_Py_C1003_b = I_NAI_H2y3z_S_C1003_b+ABY*I_NAI_Gy3z_S_C1003_b;
  Double I_NAI_G4z_Py_C1003_b = I_NAI_Hy4z_S_C1003_b+ABY*I_NAI_G4z_S_C1003_b;
  Double I_NAI_G3xz_Pz_C1003_b = I_NAI_H3x2z_S_C1003_b+ABZ*I_NAI_G3xz_S_C1003_b;
  Double I_NAI_G2xyz_Pz_C1003_b = I_NAI_H2xy2z_S_C1003_b+ABZ*I_NAI_G2xyz_S_C1003_b;
  Double I_NAI_G2x2z_Pz_C1003_b = I_NAI_H2x3z_S_C1003_b+ABZ*I_NAI_G2x2z_S_C1003_b;
  Double I_NAI_Gx2yz_Pz_C1003_b = I_NAI_Hx2y2z_S_C1003_b+ABZ*I_NAI_Gx2yz_S_C1003_b;
  Double I_NAI_Gxy2z_Pz_C1003_b = I_NAI_Hxy3z_S_C1003_b+ABZ*I_NAI_Gxy2z_S_C1003_b;
  Double I_NAI_Gx3z_Pz_C1003_b = I_NAI_Hx4z_S_C1003_b+ABZ*I_NAI_Gx3z_S_C1003_b;
  Double I_NAI_G3yz_Pz_C1003_b = I_NAI_H3y2z_S_C1003_b+ABZ*I_NAI_G3yz_S_C1003_b;
  Double I_NAI_G2y2z_Pz_C1003_b = I_NAI_H2y3z_S_C1003_b+ABZ*I_NAI_G2y2z_S_C1003_b;
  Double I_NAI_Gy3z_Pz_C1003_b = I_NAI_Hy4z_S_C1003_b+ABZ*I_NAI_Gy3z_S_C1003_b;
  Double I_NAI_G4z_Pz_C1003_b = I_NAI_H5z_S_C1003_b+ABZ*I_NAI_G4z_S_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_D_C1003_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C1003_b
   * RHS shell quartet name: SQ_NAI_F_P_C1003_b
   ************************************************************/
  Double I_NAI_F3x_D2x_C1003_b = I_NAI_G4x_Px_C1003_b+ABX*I_NAI_F3x_Px_C1003_b;
  Double I_NAI_F2xy_D2x_C1003_b = I_NAI_G3xy_Px_C1003_b+ABX*I_NAI_F2xy_Px_C1003_b;
  Double I_NAI_F2xz_D2x_C1003_b = I_NAI_G3xz_Px_C1003_b+ABX*I_NAI_F2xz_Px_C1003_b;
  Double I_NAI_Fx2y_D2x_C1003_b = I_NAI_G2x2y_Px_C1003_b+ABX*I_NAI_Fx2y_Px_C1003_b;
  Double I_NAI_Fxyz_D2x_C1003_b = I_NAI_G2xyz_Px_C1003_b+ABX*I_NAI_Fxyz_Px_C1003_b;
  Double I_NAI_Fx2z_D2x_C1003_b = I_NAI_G2x2z_Px_C1003_b+ABX*I_NAI_Fx2z_Px_C1003_b;
  Double I_NAI_F3y_D2x_C1003_b = I_NAI_Gx3y_Px_C1003_b+ABX*I_NAI_F3y_Px_C1003_b;
  Double I_NAI_F2yz_D2x_C1003_b = I_NAI_Gx2yz_Px_C1003_b+ABX*I_NAI_F2yz_Px_C1003_b;
  Double I_NAI_Fy2z_D2x_C1003_b = I_NAI_Gxy2z_Px_C1003_b+ABX*I_NAI_Fy2z_Px_C1003_b;
  Double I_NAI_F3z_D2x_C1003_b = I_NAI_Gx3z_Px_C1003_b+ABX*I_NAI_F3z_Px_C1003_b;
  Double I_NAI_F3x_Dxy_C1003_b = I_NAI_G3xy_Px_C1003_b+ABY*I_NAI_F3x_Px_C1003_b;
  Double I_NAI_F2xy_Dxy_C1003_b = I_NAI_G2x2y_Px_C1003_b+ABY*I_NAI_F2xy_Px_C1003_b;
  Double I_NAI_F2xz_Dxy_C1003_b = I_NAI_G2xyz_Px_C1003_b+ABY*I_NAI_F2xz_Px_C1003_b;
  Double I_NAI_Fx2y_Dxy_C1003_b = I_NAI_Gx3y_Px_C1003_b+ABY*I_NAI_Fx2y_Px_C1003_b;
  Double I_NAI_Fxyz_Dxy_C1003_b = I_NAI_Gx2yz_Px_C1003_b+ABY*I_NAI_Fxyz_Px_C1003_b;
  Double I_NAI_Fx2z_Dxy_C1003_b = I_NAI_Gxy2z_Px_C1003_b+ABY*I_NAI_Fx2z_Px_C1003_b;
  Double I_NAI_F3y_Dxy_C1003_b = I_NAI_G4y_Px_C1003_b+ABY*I_NAI_F3y_Px_C1003_b;
  Double I_NAI_F2yz_Dxy_C1003_b = I_NAI_G3yz_Px_C1003_b+ABY*I_NAI_F2yz_Px_C1003_b;
  Double I_NAI_Fy2z_Dxy_C1003_b = I_NAI_G2y2z_Px_C1003_b+ABY*I_NAI_Fy2z_Px_C1003_b;
  Double I_NAI_F3z_Dxy_C1003_b = I_NAI_Gy3z_Px_C1003_b+ABY*I_NAI_F3z_Px_C1003_b;
  Double I_NAI_F3x_Dxz_C1003_b = I_NAI_G3xz_Px_C1003_b+ABZ*I_NAI_F3x_Px_C1003_b;
  Double I_NAI_F2xy_Dxz_C1003_b = I_NAI_G2xyz_Px_C1003_b+ABZ*I_NAI_F2xy_Px_C1003_b;
  Double I_NAI_F2xz_Dxz_C1003_b = I_NAI_G2x2z_Px_C1003_b+ABZ*I_NAI_F2xz_Px_C1003_b;
  Double I_NAI_Fx2y_Dxz_C1003_b = I_NAI_Gx2yz_Px_C1003_b+ABZ*I_NAI_Fx2y_Px_C1003_b;
  Double I_NAI_Fxyz_Dxz_C1003_b = I_NAI_Gxy2z_Px_C1003_b+ABZ*I_NAI_Fxyz_Px_C1003_b;
  Double I_NAI_Fx2z_Dxz_C1003_b = I_NAI_Gx3z_Px_C1003_b+ABZ*I_NAI_Fx2z_Px_C1003_b;
  Double I_NAI_F3y_Dxz_C1003_b = I_NAI_G3yz_Px_C1003_b+ABZ*I_NAI_F3y_Px_C1003_b;
  Double I_NAI_F2yz_Dxz_C1003_b = I_NAI_G2y2z_Px_C1003_b+ABZ*I_NAI_F2yz_Px_C1003_b;
  Double I_NAI_Fy2z_Dxz_C1003_b = I_NAI_Gy3z_Px_C1003_b+ABZ*I_NAI_Fy2z_Px_C1003_b;
  Double I_NAI_F3z_Dxz_C1003_b = I_NAI_G4z_Px_C1003_b+ABZ*I_NAI_F3z_Px_C1003_b;
  Double I_NAI_F3x_D2y_C1003_b = I_NAI_G3xy_Py_C1003_b+ABY*I_NAI_F3x_Py_C1003_b;
  Double I_NAI_F2xy_D2y_C1003_b = I_NAI_G2x2y_Py_C1003_b+ABY*I_NAI_F2xy_Py_C1003_b;
  Double I_NAI_F2xz_D2y_C1003_b = I_NAI_G2xyz_Py_C1003_b+ABY*I_NAI_F2xz_Py_C1003_b;
  Double I_NAI_Fx2y_D2y_C1003_b = I_NAI_Gx3y_Py_C1003_b+ABY*I_NAI_Fx2y_Py_C1003_b;
  Double I_NAI_Fxyz_D2y_C1003_b = I_NAI_Gx2yz_Py_C1003_b+ABY*I_NAI_Fxyz_Py_C1003_b;
  Double I_NAI_Fx2z_D2y_C1003_b = I_NAI_Gxy2z_Py_C1003_b+ABY*I_NAI_Fx2z_Py_C1003_b;
  Double I_NAI_F3y_D2y_C1003_b = I_NAI_G4y_Py_C1003_b+ABY*I_NAI_F3y_Py_C1003_b;
  Double I_NAI_F2yz_D2y_C1003_b = I_NAI_G3yz_Py_C1003_b+ABY*I_NAI_F2yz_Py_C1003_b;
  Double I_NAI_Fy2z_D2y_C1003_b = I_NAI_G2y2z_Py_C1003_b+ABY*I_NAI_Fy2z_Py_C1003_b;
  Double I_NAI_F3z_D2y_C1003_b = I_NAI_Gy3z_Py_C1003_b+ABY*I_NAI_F3z_Py_C1003_b;
  Double I_NAI_F3x_Dyz_C1003_b = I_NAI_G3xz_Py_C1003_b+ABZ*I_NAI_F3x_Py_C1003_b;
  Double I_NAI_F2xy_Dyz_C1003_b = I_NAI_G2xyz_Py_C1003_b+ABZ*I_NAI_F2xy_Py_C1003_b;
  Double I_NAI_F2xz_Dyz_C1003_b = I_NAI_G2x2z_Py_C1003_b+ABZ*I_NAI_F2xz_Py_C1003_b;
  Double I_NAI_Fx2y_Dyz_C1003_b = I_NAI_Gx2yz_Py_C1003_b+ABZ*I_NAI_Fx2y_Py_C1003_b;
  Double I_NAI_Fxyz_Dyz_C1003_b = I_NAI_Gxy2z_Py_C1003_b+ABZ*I_NAI_Fxyz_Py_C1003_b;
  Double I_NAI_Fx2z_Dyz_C1003_b = I_NAI_Gx3z_Py_C1003_b+ABZ*I_NAI_Fx2z_Py_C1003_b;
  Double I_NAI_F3y_Dyz_C1003_b = I_NAI_G3yz_Py_C1003_b+ABZ*I_NAI_F3y_Py_C1003_b;
  Double I_NAI_F2yz_Dyz_C1003_b = I_NAI_G2y2z_Py_C1003_b+ABZ*I_NAI_F2yz_Py_C1003_b;
  Double I_NAI_Fy2z_Dyz_C1003_b = I_NAI_Gy3z_Py_C1003_b+ABZ*I_NAI_Fy2z_Py_C1003_b;
  Double I_NAI_F3z_Dyz_C1003_b = I_NAI_G4z_Py_C1003_b+ABZ*I_NAI_F3z_Py_C1003_b;
  Double I_NAI_F3x_D2z_C1003_b = I_NAI_G3xz_Pz_C1003_b+ABZ*I_NAI_F3x_Pz_C1003_b;
  Double I_NAI_F2xy_D2z_C1003_b = I_NAI_G2xyz_Pz_C1003_b+ABZ*I_NAI_F2xy_Pz_C1003_b;
  Double I_NAI_F2xz_D2z_C1003_b = I_NAI_G2x2z_Pz_C1003_b+ABZ*I_NAI_F2xz_Pz_C1003_b;
  Double I_NAI_Fx2y_D2z_C1003_b = I_NAI_Gx2yz_Pz_C1003_b+ABZ*I_NAI_Fx2y_Pz_C1003_b;
  Double I_NAI_Fxyz_D2z_C1003_b = I_NAI_Gxy2z_Pz_C1003_b+ABZ*I_NAI_Fxyz_Pz_C1003_b;
  Double I_NAI_Fx2z_D2z_C1003_b = I_NAI_Gx3z_Pz_C1003_b+ABZ*I_NAI_Fx2z_Pz_C1003_b;
  Double I_NAI_F3y_D2z_C1003_b = I_NAI_G3yz_Pz_C1003_b+ABZ*I_NAI_F3y_Pz_C1003_b;
  Double I_NAI_F2yz_D2z_C1003_b = I_NAI_G2y2z_Pz_C1003_b+ABZ*I_NAI_F2yz_Pz_C1003_b;
  Double I_NAI_Fy2z_D2z_C1003_b = I_NAI_Gy3z_Pz_C1003_b+ABZ*I_NAI_Fy2z_Pz_C1003_b;
  Double I_NAI_F3z_D2z_C1003_b = I_NAI_G4z_Pz_C1003_b+ABZ*I_NAI_F3z_Pz_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C3_a
   * RHS shell quartet name: SQ_NAI_D_S_C3
   ************************************************************/
  abcd[0] = 2.0E0*I_NAI_G4x_S_C3_a-3*I_NAI_D2x_S_C3;
  abcd[1] = 2.0E0*I_NAI_G3xy_S_C3_a-2*I_NAI_Dxy_S_C3;
  abcd[2] = 2.0E0*I_NAI_G3xz_S_C3_a-2*I_NAI_Dxz_S_C3;
  abcd[3] = 2.0E0*I_NAI_G2x2y_S_C3_a-1*I_NAI_D2y_S_C3;
  abcd[4] = 2.0E0*I_NAI_G2xyz_S_C3_a-1*I_NAI_Dyz_S_C3;
  abcd[5] = 2.0E0*I_NAI_G2x2z_S_C3_a-1*I_NAI_D2z_S_C3;
  abcd[6] = 2.0E0*I_NAI_Gx3y_S_C3_a;
  abcd[7] = 2.0E0*I_NAI_Gx2yz_S_C3_a;
  abcd[8] = 2.0E0*I_NAI_Gxy2z_S_C3_a;
  abcd[9] = 2.0E0*I_NAI_Gx3z_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C1003_a
   * RHS shell quartet name: SQ_NAI_D_P_C1003
   ************************************************************/
  abcd[10] = 2.0E0*I_NAI_G4x_Px_C1003_a-3*I_NAI_D2x_Px_C1003;
  abcd[11] = 2.0E0*I_NAI_G3xy_Px_C1003_a-2*I_NAI_Dxy_Px_C1003;
  abcd[12] = 2.0E0*I_NAI_G3xz_Px_C1003_a-2*I_NAI_Dxz_Px_C1003;
  abcd[13] = 2.0E0*I_NAI_G2x2y_Px_C1003_a-1*I_NAI_D2y_Px_C1003;
  abcd[14] = 2.0E0*I_NAI_G2xyz_Px_C1003_a-1*I_NAI_Dyz_Px_C1003;
  abcd[15] = 2.0E0*I_NAI_G2x2z_Px_C1003_a-1*I_NAI_D2z_Px_C1003;
  abcd[16] = 2.0E0*I_NAI_Gx3y_Px_C1003_a;
  abcd[17] = 2.0E0*I_NAI_Gx2yz_Px_C1003_a;
  abcd[18] = 2.0E0*I_NAI_Gxy2z_Px_C1003_a;
  abcd[19] = 2.0E0*I_NAI_Gx3z_Px_C1003_a;
  abcd[20] = 2.0E0*I_NAI_G4x_Py_C1003_a-3*I_NAI_D2x_Py_C1003;
  abcd[21] = 2.0E0*I_NAI_G3xy_Py_C1003_a-2*I_NAI_Dxy_Py_C1003;
  abcd[22] = 2.0E0*I_NAI_G3xz_Py_C1003_a-2*I_NAI_Dxz_Py_C1003;
  abcd[23] = 2.0E0*I_NAI_G2x2y_Py_C1003_a-1*I_NAI_D2y_Py_C1003;
  abcd[24] = 2.0E0*I_NAI_G2xyz_Py_C1003_a-1*I_NAI_Dyz_Py_C1003;
  abcd[25] = 2.0E0*I_NAI_G2x2z_Py_C1003_a-1*I_NAI_D2z_Py_C1003;
  abcd[26] = 2.0E0*I_NAI_Gx3y_Py_C1003_a;
  abcd[27] = 2.0E0*I_NAI_Gx2yz_Py_C1003_a;
  abcd[28] = 2.0E0*I_NAI_Gxy2z_Py_C1003_a;
  abcd[29] = 2.0E0*I_NAI_Gx3z_Py_C1003_a;
  abcd[30] = 2.0E0*I_NAI_G4x_Pz_C1003_a-3*I_NAI_D2x_Pz_C1003;
  abcd[31] = 2.0E0*I_NAI_G3xy_Pz_C1003_a-2*I_NAI_Dxy_Pz_C1003;
  abcd[32] = 2.0E0*I_NAI_G3xz_Pz_C1003_a-2*I_NAI_Dxz_Pz_C1003;
  abcd[33] = 2.0E0*I_NAI_G2x2y_Pz_C1003_a-1*I_NAI_D2y_Pz_C1003;
  abcd[34] = 2.0E0*I_NAI_G2xyz_Pz_C1003_a-1*I_NAI_Dyz_Pz_C1003;
  abcd[35] = 2.0E0*I_NAI_G2x2z_Pz_C1003_a-1*I_NAI_D2z_Pz_C1003;
  abcd[36] = 2.0E0*I_NAI_Gx3y_Pz_C1003_a;
  abcd[37] = 2.0E0*I_NAI_Gx2yz_Pz_C1003_a;
  abcd[38] = 2.0E0*I_NAI_Gxy2z_Pz_C1003_a;
  abcd[39] = 2.0E0*I_NAI_Gx3z_Pz_C1003_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C3_a
   * RHS shell quartet name: SQ_NAI_D_S_C3
   ************************************************************/
  abcd[40] = 2.0E0*I_NAI_G3xy_S_C3_a;
  abcd[41] = 2.0E0*I_NAI_G2x2y_S_C3_a-1*I_NAI_D2x_S_C3;
  abcd[42] = 2.0E0*I_NAI_G2xyz_S_C3_a;
  abcd[43] = 2.0E0*I_NAI_Gx3y_S_C3_a-2*I_NAI_Dxy_S_C3;
  abcd[44] = 2.0E0*I_NAI_Gx2yz_S_C3_a-1*I_NAI_Dxz_S_C3;
  abcd[45] = 2.0E0*I_NAI_Gxy2z_S_C3_a;
  abcd[46] = 2.0E0*I_NAI_G4y_S_C3_a-3*I_NAI_D2y_S_C3;
  abcd[47] = 2.0E0*I_NAI_G3yz_S_C3_a-2*I_NAI_Dyz_S_C3;
  abcd[48] = 2.0E0*I_NAI_G2y2z_S_C3_a-1*I_NAI_D2z_S_C3;
  abcd[49] = 2.0E0*I_NAI_Gy3z_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C1003_a
   * RHS shell quartet name: SQ_NAI_D_P_C1003
   ************************************************************/
  abcd[50] = 2.0E0*I_NAI_G3xy_Px_C1003_a;
  abcd[51] = 2.0E0*I_NAI_G2x2y_Px_C1003_a-1*I_NAI_D2x_Px_C1003;
  abcd[52] = 2.0E0*I_NAI_G2xyz_Px_C1003_a;
  abcd[53] = 2.0E0*I_NAI_Gx3y_Px_C1003_a-2*I_NAI_Dxy_Px_C1003;
  abcd[54] = 2.0E0*I_NAI_Gx2yz_Px_C1003_a-1*I_NAI_Dxz_Px_C1003;
  abcd[55] = 2.0E0*I_NAI_Gxy2z_Px_C1003_a;
  abcd[56] = 2.0E0*I_NAI_G4y_Px_C1003_a-3*I_NAI_D2y_Px_C1003;
  abcd[57] = 2.0E0*I_NAI_G3yz_Px_C1003_a-2*I_NAI_Dyz_Px_C1003;
  abcd[58] = 2.0E0*I_NAI_G2y2z_Px_C1003_a-1*I_NAI_D2z_Px_C1003;
  abcd[59] = 2.0E0*I_NAI_Gy3z_Px_C1003_a;
  abcd[60] = 2.0E0*I_NAI_G3xy_Py_C1003_a;
  abcd[61] = 2.0E0*I_NAI_G2x2y_Py_C1003_a-1*I_NAI_D2x_Py_C1003;
  abcd[62] = 2.0E0*I_NAI_G2xyz_Py_C1003_a;
  abcd[63] = 2.0E0*I_NAI_Gx3y_Py_C1003_a-2*I_NAI_Dxy_Py_C1003;
  abcd[64] = 2.0E0*I_NAI_Gx2yz_Py_C1003_a-1*I_NAI_Dxz_Py_C1003;
  abcd[65] = 2.0E0*I_NAI_Gxy2z_Py_C1003_a;
  abcd[66] = 2.0E0*I_NAI_G4y_Py_C1003_a-3*I_NAI_D2y_Py_C1003;
  abcd[67] = 2.0E0*I_NAI_G3yz_Py_C1003_a-2*I_NAI_Dyz_Py_C1003;
  abcd[68] = 2.0E0*I_NAI_G2y2z_Py_C1003_a-1*I_NAI_D2z_Py_C1003;
  abcd[69] = 2.0E0*I_NAI_Gy3z_Py_C1003_a;
  abcd[70] = 2.0E0*I_NAI_G3xy_Pz_C1003_a;
  abcd[71] = 2.0E0*I_NAI_G2x2y_Pz_C1003_a-1*I_NAI_D2x_Pz_C1003;
  abcd[72] = 2.0E0*I_NAI_G2xyz_Pz_C1003_a;
  abcd[73] = 2.0E0*I_NAI_Gx3y_Pz_C1003_a-2*I_NAI_Dxy_Pz_C1003;
  abcd[74] = 2.0E0*I_NAI_Gx2yz_Pz_C1003_a-1*I_NAI_Dxz_Pz_C1003;
  abcd[75] = 2.0E0*I_NAI_Gxy2z_Pz_C1003_a;
  abcd[76] = 2.0E0*I_NAI_G4y_Pz_C1003_a-3*I_NAI_D2y_Pz_C1003;
  abcd[77] = 2.0E0*I_NAI_G3yz_Pz_C1003_a-2*I_NAI_Dyz_Pz_C1003;
  abcd[78] = 2.0E0*I_NAI_G2y2z_Pz_C1003_a-1*I_NAI_D2z_Pz_C1003;
  abcd[79] = 2.0E0*I_NAI_Gy3z_Pz_C1003_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C3_a
   * RHS shell quartet name: SQ_NAI_D_S_C3
   ************************************************************/
  abcd[80] = 2.0E0*I_NAI_G3xz_S_C3_a;
  abcd[81] = 2.0E0*I_NAI_G2xyz_S_C3_a;
  abcd[82] = 2.0E0*I_NAI_G2x2z_S_C3_a-1*I_NAI_D2x_S_C3;
  abcd[83] = 2.0E0*I_NAI_Gx2yz_S_C3_a;
  abcd[84] = 2.0E0*I_NAI_Gxy2z_S_C3_a-1*I_NAI_Dxy_S_C3;
  abcd[85] = 2.0E0*I_NAI_Gx3z_S_C3_a-2*I_NAI_Dxz_S_C3;
  abcd[86] = 2.0E0*I_NAI_G3yz_S_C3_a;
  abcd[87] = 2.0E0*I_NAI_G2y2z_S_C3_a-1*I_NAI_D2y_S_C3;
  abcd[88] = 2.0E0*I_NAI_Gy3z_S_C3_a-2*I_NAI_Dyz_S_C3;
  abcd[89] = 2.0E0*I_NAI_G4z_S_C3_a-3*I_NAI_D2z_S_C3;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_C1003_a
   * RHS shell quartet name: SQ_NAI_D_P_C1003
   ************************************************************/
  abcd[90] = 2.0E0*I_NAI_G3xz_Px_C1003_a;
  abcd[91] = 2.0E0*I_NAI_G2xyz_Px_C1003_a;
  abcd[92] = 2.0E0*I_NAI_G2x2z_Px_C1003_a-1*I_NAI_D2x_Px_C1003;
  abcd[93] = 2.0E0*I_NAI_Gx2yz_Px_C1003_a;
  abcd[94] = 2.0E0*I_NAI_Gxy2z_Px_C1003_a-1*I_NAI_Dxy_Px_C1003;
  abcd[95] = 2.0E0*I_NAI_Gx3z_Px_C1003_a-2*I_NAI_Dxz_Px_C1003;
  abcd[96] = 2.0E0*I_NAI_G3yz_Px_C1003_a;
  abcd[97] = 2.0E0*I_NAI_G2y2z_Px_C1003_a-1*I_NAI_D2y_Px_C1003;
  abcd[98] = 2.0E0*I_NAI_Gy3z_Px_C1003_a-2*I_NAI_Dyz_Px_C1003;
  abcd[99] = 2.0E0*I_NAI_G4z_Px_C1003_a-3*I_NAI_D2z_Px_C1003;
  abcd[100] = 2.0E0*I_NAI_G3xz_Py_C1003_a;
  abcd[101] = 2.0E0*I_NAI_G2xyz_Py_C1003_a;
  abcd[102] = 2.0E0*I_NAI_G2x2z_Py_C1003_a-1*I_NAI_D2x_Py_C1003;
  abcd[103] = 2.0E0*I_NAI_Gx2yz_Py_C1003_a;
  abcd[104] = 2.0E0*I_NAI_Gxy2z_Py_C1003_a-1*I_NAI_Dxy_Py_C1003;
  abcd[105] = 2.0E0*I_NAI_Gx3z_Py_C1003_a-2*I_NAI_Dxz_Py_C1003;
  abcd[106] = 2.0E0*I_NAI_G3yz_Py_C1003_a;
  abcd[107] = 2.0E0*I_NAI_G2y2z_Py_C1003_a-1*I_NAI_D2y_Py_C1003;
  abcd[108] = 2.0E0*I_NAI_Gy3z_Py_C1003_a-2*I_NAI_Dyz_Py_C1003;
  abcd[109] = 2.0E0*I_NAI_G4z_Py_C1003_a-3*I_NAI_D2z_Py_C1003;
  abcd[110] = 2.0E0*I_NAI_G3xz_Pz_C1003_a;
  abcd[111] = 2.0E0*I_NAI_G2xyz_Pz_C1003_a;
  abcd[112] = 2.0E0*I_NAI_G2x2z_Pz_C1003_a-1*I_NAI_D2x_Pz_C1003;
  abcd[113] = 2.0E0*I_NAI_Gx2yz_Pz_C1003_a;
  abcd[114] = 2.0E0*I_NAI_Gxy2z_Pz_C1003_a-1*I_NAI_Dxy_Pz_C1003;
  abcd[115] = 2.0E0*I_NAI_Gx3z_Pz_C1003_a-2*I_NAI_Dxz_Pz_C1003;
  abcd[116] = 2.0E0*I_NAI_G3yz_Pz_C1003_a;
  abcd[117] = 2.0E0*I_NAI_G2y2z_Pz_C1003_a-1*I_NAI_D2y_Pz_C1003;
  abcd[118] = 2.0E0*I_NAI_Gy3z_Pz_C1003_a-2*I_NAI_Dyz_Pz_C1003;
  abcd[119] = 2.0E0*I_NAI_G4z_Pz_C1003_a-3*I_NAI_D2z_Pz_C1003;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C3_b
   ************************************************************/
  abcd[120] = 2.0E0*I_NAI_F3x_Px_C3_b;
  abcd[121] = 2.0E0*I_NAI_F2xy_Px_C3_b;
  abcd[122] = 2.0E0*I_NAI_F2xz_Px_C3_b;
  abcd[123] = 2.0E0*I_NAI_Fx2y_Px_C3_b;
  abcd[124] = 2.0E0*I_NAI_Fxyz_Px_C3_b;
  abcd[125] = 2.0E0*I_NAI_Fx2z_Px_C3_b;
  abcd[126] = 2.0E0*I_NAI_F3y_Px_C3_b;
  abcd[127] = 2.0E0*I_NAI_F2yz_Px_C3_b;
  abcd[128] = 2.0E0*I_NAI_Fy2z_Px_C3_b;
  abcd[129] = 2.0E0*I_NAI_F3z_Px_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D_C1003_b
   * RHS shell quartet name: SQ_NAI_F_S_C1003
   ************************************************************/
  abcd[130] = 2.0E0*I_NAI_F3x_D2x_C1003_b-1*I_NAI_F3x_S_C1003;
  abcd[131] = 2.0E0*I_NAI_F2xy_D2x_C1003_b-1*I_NAI_F2xy_S_C1003;
  abcd[132] = 2.0E0*I_NAI_F2xz_D2x_C1003_b-1*I_NAI_F2xz_S_C1003;
  abcd[133] = 2.0E0*I_NAI_Fx2y_D2x_C1003_b-1*I_NAI_Fx2y_S_C1003;
  abcd[134] = 2.0E0*I_NAI_Fxyz_D2x_C1003_b-1*I_NAI_Fxyz_S_C1003;
  abcd[135] = 2.0E0*I_NAI_Fx2z_D2x_C1003_b-1*I_NAI_Fx2z_S_C1003;
  abcd[136] = 2.0E0*I_NAI_F3y_D2x_C1003_b-1*I_NAI_F3y_S_C1003;
  abcd[137] = 2.0E0*I_NAI_F2yz_D2x_C1003_b-1*I_NAI_F2yz_S_C1003;
  abcd[138] = 2.0E0*I_NAI_Fy2z_D2x_C1003_b-1*I_NAI_Fy2z_S_C1003;
  abcd[139] = 2.0E0*I_NAI_F3z_D2x_C1003_b-1*I_NAI_F3z_S_C1003;
  abcd[140] = 2.0E0*I_NAI_F3x_Dxy_C1003_b;
  abcd[141] = 2.0E0*I_NAI_F2xy_Dxy_C1003_b;
  abcd[142] = 2.0E0*I_NAI_F2xz_Dxy_C1003_b;
  abcd[143] = 2.0E0*I_NAI_Fx2y_Dxy_C1003_b;
  abcd[144] = 2.0E0*I_NAI_Fxyz_Dxy_C1003_b;
  abcd[145] = 2.0E0*I_NAI_Fx2z_Dxy_C1003_b;
  abcd[146] = 2.0E0*I_NAI_F3y_Dxy_C1003_b;
  abcd[147] = 2.0E0*I_NAI_F2yz_Dxy_C1003_b;
  abcd[148] = 2.0E0*I_NAI_Fy2z_Dxy_C1003_b;
  abcd[149] = 2.0E0*I_NAI_F3z_Dxy_C1003_b;
  abcd[150] = 2.0E0*I_NAI_F3x_Dxz_C1003_b;
  abcd[151] = 2.0E0*I_NAI_F2xy_Dxz_C1003_b;
  abcd[152] = 2.0E0*I_NAI_F2xz_Dxz_C1003_b;
  abcd[153] = 2.0E0*I_NAI_Fx2y_Dxz_C1003_b;
  abcd[154] = 2.0E0*I_NAI_Fxyz_Dxz_C1003_b;
  abcd[155] = 2.0E0*I_NAI_Fx2z_Dxz_C1003_b;
  abcd[156] = 2.0E0*I_NAI_F3y_Dxz_C1003_b;
  abcd[157] = 2.0E0*I_NAI_F2yz_Dxz_C1003_b;
  abcd[158] = 2.0E0*I_NAI_Fy2z_Dxz_C1003_b;
  abcd[159] = 2.0E0*I_NAI_F3z_Dxz_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C3_b
   ************************************************************/
  abcd[160] = 2.0E0*I_NAI_F3x_Py_C3_b;
  abcd[161] = 2.0E0*I_NAI_F2xy_Py_C3_b;
  abcd[162] = 2.0E0*I_NAI_F2xz_Py_C3_b;
  abcd[163] = 2.0E0*I_NAI_Fx2y_Py_C3_b;
  abcd[164] = 2.0E0*I_NAI_Fxyz_Py_C3_b;
  abcd[165] = 2.0E0*I_NAI_Fx2z_Py_C3_b;
  abcd[166] = 2.0E0*I_NAI_F3y_Py_C3_b;
  abcd[167] = 2.0E0*I_NAI_F2yz_Py_C3_b;
  abcd[168] = 2.0E0*I_NAI_Fy2z_Py_C3_b;
  abcd[169] = 2.0E0*I_NAI_F3z_Py_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D_C1003_b
   * RHS shell quartet name: SQ_NAI_F_S_C1003
   ************************************************************/
  abcd[170] = 2.0E0*I_NAI_F3x_Dxy_C1003_b;
  abcd[171] = 2.0E0*I_NAI_F2xy_Dxy_C1003_b;
  abcd[172] = 2.0E0*I_NAI_F2xz_Dxy_C1003_b;
  abcd[173] = 2.0E0*I_NAI_Fx2y_Dxy_C1003_b;
  abcd[174] = 2.0E0*I_NAI_Fxyz_Dxy_C1003_b;
  abcd[175] = 2.0E0*I_NAI_Fx2z_Dxy_C1003_b;
  abcd[176] = 2.0E0*I_NAI_F3y_Dxy_C1003_b;
  abcd[177] = 2.0E0*I_NAI_F2yz_Dxy_C1003_b;
  abcd[178] = 2.0E0*I_NAI_Fy2z_Dxy_C1003_b;
  abcd[179] = 2.0E0*I_NAI_F3z_Dxy_C1003_b;
  abcd[180] = 2.0E0*I_NAI_F3x_D2y_C1003_b-1*I_NAI_F3x_S_C1003;
  abcd[181] = 2.0E0*I_NAI_F2xy_D2y_C1003_b-1*I_NAI_F2xy_S_C1003;
  abcd[182] = 2.0E0*I_NAI_F2xz_D2y_C1003_b-1*I_NAI_F2xz_S_C1003;
  abcd[183] = 2.0E0*I_NAI_Fx2y_D2y_C1003_b-1*I_NAI_Fx2y_S_C1003;
  abcd[184] = 2.0E0*I_NAI_Fxyz_D2y_C1003_b-1*I_NAI_Fxyz_S_C1003;
  abcd[185] = 2.0E0*I_NAI_Fx2z_D2y_C1003_b-1*I_NAI_Fx2z_S_C1003;
  abcd[186] = 2.0E0*I_NAI_F3y_D2y_C1003_b-1*I_NAI_F3y_S_C1003;
  abcd[187] = 2.0E0*I_NAI_F2yz_D2y_C1003_b-1*I_NAI_F2yz_S_C1003;
  abcd[188] = 2.0E0*I_NAI_Fy2z_D2y_C1003_b-1*I_NAI_Fy2z_S_C1003;
  abcd[189] = 2.0E0*I_NAI_F3z_D2y_C1003_b-1*I_NAI_F3z_S_C1003;
  abcd[190] = 2.0E0*I_NAI_F3x_Dyz_C1003_b;
  abcd[191] = 2.0E0*I_NAI_F2xy_Dyz_C1003_b;
  abcd[192] = 2.0E0*I_NAI_F2xz_Dyz_C1003_b;
  abcd[193] = 2.0E0*I_NAI_Fx2y_Dyz_C1003_b;
  abcd[194] = 2.0E0*I_NAI_Fxyz_Dyz_C1003_b;
  abcd[195] = 2.0E0*I_NAI_Fx2z_Dyz_C1003_b;
  abcd[196] = 2.0E0*I_NAI_F3y_Dyz_C1003_b;
  abcd[197] = 2.0E0*I_NAI_F2yz_Dyz_C1003_b;
  abcd[198] = 2.0E0*I_NAI_Fy2z_Dyz_C1003_b;
  abcd[199] = 2.0E0*I_NAI_F3z_Dyz_C1003_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_S_C3_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C3_b
   ************************************************************/
  abcd[200] = 2.0E0*I_NAI_F3x_Pz_C3_b;
  abcd[201] = 2.0E0*I_NAI_F2xy_Pz_C3_b;
  abcd[202] = 2.0E0*I_NAI_F2xz_Pz_C3_b;
  abcd[203] = 2.0E0*I_NAI_Fx2y_Pz_C3_b;
  abcd[204] = 2.0E0*I_NAI_Fxyz_Pz_C3_b;
  abcd[205] = 2.0E0*I_NAI_Fx2z_Pz_C3_b;
  abcd[206] = 2.0E0*I_NAI_F3y_Pz_C3_b;
  abcd[207] = 2.0E0*I_NAI_F2yz_Pz_C3_b;
  abcd[208] = 2.0E0*I_NAI_Fy2z_Pz_C3_b;
  abcd[209] = 2.0E0*I_NAI_F3z_Pz_C3_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1003_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D_C1003_b
   * RHS shell quartet name: SQ_NAI_F_S_C1003
   ************************************************************/
  abcd[210] = 2.0E0*I_NAI_F3x_Dxz_C1003_b;
  abcd[211] = 2.0E0*I_NAI_F2xy_Dxz_C1003_b;
  abcd[212] = 2.0E0*I_NAI_F2xz_Dxz_C1003_b;
  abcd[213] = 2.0E0*I_NAI_Fx2y_Dxz_C1003_b;
  abcd[214] = 2.0E0*I_NAI_Fxyz_Dxz_C1003_b;
  abcd[215] = 2.0E0*I_NAI_Fx2z_Dxz_C1003_b;
  abcd[216] = 2.0E0*I_NAI_F3y_Dxz_C1003_b;
  abcd[217] = 2.0E0*I_NAI_F2yz_Dxz_C1003_b;
  abcd[218] = 2.0E0*I_NAI_Fy2z_Dxz_C1003_b;
  abcd[219] = 2.0E0*I_NAI_F3z_Dxz_C1003_b;
  abcd[220] = 2.0E0*I_NAI_F3x_Dyz_C1003_b;
  abcd[221] = 2.0E0*I_NAI_F2xy_Dyz_C1003_b;
  abcd[222] = 2.0E0*I_NAI_F2xz_Dyz_C1003_b;
  abcd[223] = 2.0E0*I_NAI_Fx2y_Dyz_C1003_b;
  abcd[224] = 2.0E0*I_NAI_Fxyz_Dyz_C1003_b;
  abcd[225] = 2.0E0*I_NAI_Fx2z_Dyz_C1003_b;
  abcd[226] = 2.0E0*I_NAI_F3y_Dyz_C1003_b;
  abcd[227] = 2.0E0*I_NAI_F2yz_Dyz_C1003_b;
  abcd[228] = 2.0E0*I_NAI_Fy2z_Dyz_C1003_b;
  abcd[229] = 2.0E0*I_NAI_F3z_Dyz_C1003_b;
  abcd[230] = 2.0E0*I_NAI_F3x_D2z_C1003_b-1*I_NAI_F3x_S_C1003;
  abcd[231] = 2.0E0*I_NAI_F2xy_D2z_C1003_b-1*I_NAI_F2xy_S_C1003;
  abcd[232] = 2.0E0*I_NAI_F2xz_D2z_C1003_b-1*I_NAI_F2xz_S_C1003;
  abcd[233] = 2.0E0*I_NAI_Fx2y_D2z_C1003_b-1*I_NAI_Fx2y_S_C1003;
  abcd[234] = 2.0E0*I_NAI_Fxyz_D2z_C1003_b-1*I_NAI_Fxyz_S_C1003;
  abcd[235] = 2.0E0*I_NAI_Fx2z_D2z_C1003_b-1*I_NAI_Fx2z_S_C1003;
  abcd[236] = 2.0E0*I_NAI_F3y_D2z_C1003_b-1*I_NAI_F3y_S_C1003;
  abcd[237] = 2.0E0*I_NAI_F2yz_D2z_C1003_b-1*I_NAI_F2yz_S_C1003;
  abcd[238] = 2.0E0*I_NAI_Fy2z_D2z_C1003_b-1*I_NAI_Fy2z_S_C1003;
  abcd[239] = 2.0E0*I_NAI_F3z_D2z_C1003_b-1*I_NAI_F3z_S_C1003;
}
