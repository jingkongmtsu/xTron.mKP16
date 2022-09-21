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

void hgp_os_nai_d_sp_d1(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_NAI_F3x_S_C2_a = 0.0E0;
  Double I_NAI_F2xy_S_C2_a = 0.0E0;
  Double I_NAI_F2xz_S_C2_a = 0.0E0;
  Double I_NAI_Fx2y_S_C2_a = 0.0E0;
  Double I_NAI_Fxyz_S_C2_a = 0.0E0;
  Double I_NAI_Fx2z_S_C2_a = 0.0E0;
  Double I_NAI_F3y_S_C2_a = 0.0E0;
  Double I_NAI_F2yz_S_C2_a = 0.0E0;
  Double I_NAI_Fy2z_S_C2_a = 0.0E0;
  Double I_NAI_F3z_S_C2_a = 0.0E0;
  Double I_NAI_Px_S_C2 = 0.0E0;
  Double I_NAI_Py_S_C2 = 0.0E0;
  Double I_NAI_Pz_S_C2 = 0.0E0;
  Double I_NAI_D2x_S_C1002 = 0.0E0;
  Double I_NAI_Dxy_S_C1002 = 0.0E0;
  Double I_NAI_Dxz_S_C1002 = 0.0E0;
  Double I_NAI_D2y_S_C1002 = 0.0E0;
  Double I_NAI_Dyz_S_C1002 = 0.0E0;
  Double I_NAI_D2z_S_C1002 = 0.0E0;
  Double I_NAI_G4x_S_C1002_a = 0.0E0;
  Double I_NAI_G3xy_S_C1002_a = 0.0E0;
  Double I_NAI_G3xz_S_C1002_a = 0.0E0;
  Double I_NAI_G2x2y_S_C1002_a = 0.0E0;
  Double I_NAI_G2xyz_S_C1002_a = 0.0E0;
  Double I_NAI_G2x2z_S_C1002_a = 0.0E0;
  Double I_NAI_Gx3y_S_C1002_a = 0.0E0;
  Double I_NAI_Gx2yz_S_C1002_a = 0.0E0;
  Double I_NAI_Gxy2z_S_C1002_a = 0.0E0;
  Double I_NAI_Gx3z_S_C1002_a = 0.0E0;
  Double I_NAI_G4y_S_C1002_a = 0.0E0;
  Double I_NAI_G3yz_S_C1002_a = 0.0E0;
  Double I_NAI_G2y2z_S_C1002_a = 0.0E0;
  Double I_NAI_Gy3z_S_C1002_a = 0.0E0;
  Double I_NAI_G4z_S_C1002_a = 0.0E0;
  Double I_NAI_F3x_S_C1002_a = 0.0E0;
  Double I_NAI_F2xy_S_C1002_a = 0.0E0;
  Double I_NAI_F2xz_S_C1002_a = 0.0E0;
  Double I_NAI_Fx2y_S_C1002_a = 0.0E0;
  Double I_NAI_Fxyz_S_C1002_a = 0.0E0;
  Double I_NAI_Fx2z_S_C1002_a = 0.0E0;
  Double I_NAI_F3y_S_C1002_a = 0.0E0;
  Double I_NAI_F2yz_S_C1002_a = 0.0E0;
  Double I_NAI_Fy2z_S_C1002_a = 0.0E0;
  Double I_NAI_F3z_S_C1002_a = 0.0E0;
  Double I_NAI_Px_S_C1002 = 0.0E0;
  Double I_NAI_Py_S_C1002 = 0.0E0;
  Double I_NAI_Pz_S_C1002 = 0.0E0;
  Double I_NAI_F3x_S_C2_b = 0.0E0;
  Double I_NAI_F2xy_S_C2_b = 0.0E0;
  Double I_NAI_F2xz_S_C2_b = 0.0E0;
  Double I_NAI_Fx2y_S_C2_b = 0.0E0;
  Double I_NAI_Fxyz_S_C2_b = 0.0E0;
  Double I_NAI_Fx2z_S_C2_b = 0.0E0;
  Double I_NAI_F3y_S_C2_b = 0.0E0;
  Double I_NAI_F2yz_S_C2_b = 0.0E0;
  Double I_NAI_Fy2z_S_C2_b = 0.0E0;
  Double I_NAI_F3z_S_C2_b = 0.0E0;
  Double I_NAI_D2x_S_C2_b = 0.0E0;
  Double I_NAI_Dxy_S_C2_b = 0.0E0;
  Double I_NAI_Dxz_S_C2_b = 0.0E0;
  Double I_NAI_D2y_S_C2_b = 0.0E0;
  Double I_NAI_Dyz_S_C2_b = 0.0E0;
  Double I_NAI_D2z_S_C2_b = 0.0E0;
  Double I_NAI_G4x_S_C1002_b = 0.0E0;
  Double I_NAI_G3xy_S_C1002_b = 0.0E0;
  Double I_NAI_G3xz_S_C1002_b = 0.0E0;
  Double I_NAI_G2x2y_S_C1002_b = 0.0E0;
  Double I_NAI_G2xyz_S_C1002_b = 0.0E0;
  Double I_NAI_G2x2z_S_C1002_b = 0.0E0;
  Double I_NAI_Gx3y_S_C1002_b = 0.0E0;
  Double I_NAI_Gx2yz_S_C1002_b = 0.0E0;
  Double I_NAI_Gxy2z_S_C1002_b = 0.0E0;
  Double I_NAI_Gx3z_S_C1002_b = 0.0E0;
  Double I_NAI_G4y_S_C1002_b = 0.0E0;
  Double I_NAI_G3yz_S_C1002_b = 0.0E0;
  Double I_NAI_G2y2z_S_C1002_b = 0.0E0;
  Double I_NAI_Gy3z_S_C1002_b = 0.0E0;
  Double I_NAI_G4z_S_C1002_b = 0.0E0;
  Double I_NAI_F3x_S_C1002_b = 0.0E0;
  Double I_NAI_F2xy_S_C1002_b = 0.0E0;
  Double I_NAI_F2xz_S_C1002_b = 0.0E0;
  Double I_NAI_Fx2y_S_C1002_b = 0.0E0;
  Double I_NAI_Fxyz_S_C1002_b = 0.0E0;
  Double I_NAI_Fx2z_S_C1002_b = 0.0E0;
  Double I_NAI_F3y_S_C1002_b = 0.0E0;
  Double I_NAI_F2yz_S_C1002_b = 0.0E0;
  Double I_NAI_Fy2z_S_C1002_b = 0.0E0;
  Double I_NAI_F3z_S_C1002_b = 0.0E0;
  Double I_NAI_D2x_S_C1002_b = 0.0E0;
  Double I_NAI_Dxy_S_C1002_b = 0.0E0;
  Double I_NAI_Dxz_S_C1002_b = 0.0E0;
  Double I_NAI_D2y_S_C1002_b = 0.0E0;
  Double I_NAI_Dyz_S_C1002_b = 0.0E0;
  Double I_NAI_D2z_S_C1002_b = 0.0E0;

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

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_NAI_S_S_vrr_d  = 0.0E0;
      double I_NAI_S_S_M1_vrr_d  = 0.0E0;
      double I_NAI_S_S_M2_vrr_d  = 0.0E0;
      double I_NAI_S_S_M3_vrr_d  = 0.0E0;
      double I_NAI_S_S_M4_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER43;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER41*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER39*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER37*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER35*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER33*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER31*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER29*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER27*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER25*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER23*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER21*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER19*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER17*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER15*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER13*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = 1.0E0+u2*ONEOVER11*I_NAI_S_S_M4_vrr;
        I_NAI_S_S_M4_vrr = ONEOVER9*I_NAI_S_S_M4_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M4_vrr  = f*I_NAI_S_S_M4_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
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

        // write the double result back to the float var
        I_NAI_S_S_vrr = static_cast<Double>(I_NAI_S_S_vrr_d);
        I_NAI_S_S_M1_vrr = static_cast<Double>(I_NAI_S_S_M1_vrr_d);
        I_NAI_S_S_M2_vrr = static_cast<Double>(I_NAI_S_S_M2_vrr_d);
        I_NAI_S_S_M3_vrr = static_cast<Double>(I_NAI_S_S_M3_vrr_d);
        I_NAI_S_S_M4_vrr = static_cast<Double>(I_NAI_S_S_M4_vrr_d);

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

#endif

      }


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
       * shell quartet name: SQ_NAI_F_S_C2_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C2_a_coefs = ic2*alpha;
      I_NAI_F3x_S_C2_a += SQ_NAI_F_S_C2_a_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C2_a += SQ_NAI_F_S_C2_a_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C2_a += SQ_NAI_F_S_C2_a_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C2_a += SQ_NAI_F_S_C2_a_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C2_a += SQ_NAI_F_S_C2_a_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C2_a += SQ_NAI_F_S_C2_a_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C2_a += SQ_NAI_F_S_C2_a_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C2_a += SQ_NAI_F_S_C2_a_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C2_a += SQ_NAI_F_S_C2_a_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C2_a += SQ_NAI_F_S_C2_a_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C2
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C2_coefs = ic2;
      I_NAI_Px_S_C2 += SQ_NAI_P_S_C2_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C2 += SQ_NAI_P_S_C2_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C2 += SQ_NAI_P_S_C2_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1002_coefs = ic2_1;
      I_NAI_D2x_S_C1002 += SQ_NAI_D_S_C1002_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1002 += SQ_NAI_D_S_C1002_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1002 += SQ_NAI_D_S_C1002_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1002 += SQ_NAI_D_S_C1002_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1002 += SQ_NAI_D_S_C1002_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1002 += SQ_NAI_D_S_C1002_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1002_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1002_a_coefs = ic2_1*alpha;
      I_NAI_G4x_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1002_a += SQ_NAI_G_S_C1002_a_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1002_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1002_a_coefs = ic2_1*alpha;
      I_NAI_F3x_S_C1002_a += SQ_NAI_F_S_C1002_a_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1002_a += SQ_NAI_F_S_C1002_a_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1002_a += SQ_NAI_F_S_C1002_a_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1002_a += SQ_NAI_F_S_C1002_a_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1002_a += SQ_NAI_F_S_C1002_a_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1002_a += SQ_NAI_F_S_C1002_a_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1002_a += SQ_NAI_F_S_C1002_a_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1002_a += SQ_NAI_F_S_C1002_a_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1002_a += SQ_NAI_F_S_C1002_a_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1002_a += SQ_NAI_F_S_C1002_a_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_C1002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_P_S_C1002_coefs = ic2_1;
      I_NAI_Px_S_C1002 += SQ_NAI_P_S_C1002_coefs*I_NAI_Px_S_vrr;
      I_NAI_Py_S_C1002 += SQ_NAI_P_S_C1002_coefs*I_NAI_Py_S_vrr;
      I_NAI_Pz_S_C1002 += SQ_NAI_P_S_C1002_coefs*I_NAI_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C2_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C2_b_coefs = ic2*beta;
      I_NAI_F3x_S_C2_b += SQ_NAI_F_S_C2_b_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C2_b += SQ_NAI_F_S_C2_b_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C2_b += SQ_NAI_F_S_C2_b_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C2_b += SQ_NAI_F_S_C2_b_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C2_b += SQ_NAI_F_S_C2_b_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C2_b += SQ_NAI_F_S_C2_b_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C2_b += SQ_NAI_F_S_C2_b_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C2_b += SQ_NAI_F_S_C2_b_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C2_b += SQ_NAI_F_S_C2_b_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C2_b += SQ_NAI_F_S_C2_b_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C2_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C2_b_coefs = ic2*beta;
      I_NAI_D2x_S_C2_b += SQ_NAI_D_S_C2_b_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C2_b += SQ_NAI_D_S_C2_b_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C2_b += SQ_NAI_D_S_C2_b_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C2_b += SQ_NAI_D_S_C2_b_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C2_b += SQ_NAI_D_S_C2_b_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C2_b += SQ_NAI_D_S_C2_b_coefs*I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_C1002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_C1002_b_coefs = ic2_1*beta;
      I_NAI_G4x_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_C1002_b += SQ_NAI_G_S_C1002_b_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_C1002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_C1002_b_coefs = ic2_1*beta;
      I_NAI_F3x_S_C1002_b += SQ_NAI_F_S_C1002_b_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_C1002_b += SQ_NAI_F_S_C1002_b_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_C1002_b += SQ_NAI_F_S_C1002_b_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_C1002_b += SQ_NAI_F_S_C1002_b_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_C1002_b += SQ_NAI_F_S_C1002_b_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_C1002_b += SQ_NAI_F_S_C1002_b_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_C1002_b += SQ_NAI_F_S_C1002_b_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_C1002_b += SQ_NAI_F_S_C1002_b_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_C1002_b += SQ_NAI_F_S_C1002_b_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_C1002_b += SQ_NAI_F_S_C1002_b_coefs*I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_C1002_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_D_S_C1002_b_coefs = ic2_1*beta;
      I_NAI_D2x_S_C1002_b += SQ_NAI_D_S_C1002_b_coefs*I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S_C1002_b += SQ_NAI_D_S_C1002_b_coefs*I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S_C1002_b += SQ_NAI_D_S_C1002_b_coefs*I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S_C1002_b += SQ_NAI_D_S_C1002_b_coefs*I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S_C1002_b += SQ_NAI_D_S_C1002_b_coefs*I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S_C1002_b += SQ_NAI_D_S_C1002_b_coefs*I_NAI_D2z_S_vrr;
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
   * shell quartet name: SQ_NAI_P_P_C1002
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_S_C1002
   * RHS shell quartet name: SQ_NAI_P_S_C1002
   ************************************************************/
  Double I_NAI_Px_Px_C1002 = I_NAI_D2x_S_C1002+ABX*I_NAI_Px_S_C1002;
  Double I_NAI_Py_Px_C1002 = I_NAI_Dxy_S_C1002+ABX*I_NAI_Py_S_C1002;
  Double I_NAI_Pz_Px_C1002 = I_NAI_Dxz_S_C1002+ABX*I_NAI_Pz_S_C1002;
  Double I_NAI_Px_Py_C1002 = I_NAI_Dxy_S_C1002+ABY*I_NAI_Px_S_C1002;
  Double I_NAI_Py_Py_C1002 = I_NAI_D2y_S_C1002+ABY*I_NAI_Py_S_C1002;
  Double I_NAI_Pz_Py_C1002 = I_NAI_Dyz_S_C1002+ABY*I_NAI_Pz_S_C1002;
  Double I_NAI_Px_Pz_C1002 = I_NAI_Dxz_S_C1002+ABZ*I_NAI_Px_S_C1002;
  Double I_NAI_Py_Pz_C1002 = I_NAI_Dyz_S_C1002+ABZ*I_NAI_Py_S_C1002;
  Double I_NAI_Pz_Pz_C1002 = I_NAI_D2z_S_C1002+ABZ*I_NAI_Pz_S_C1002;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1002_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C1002_a
   * RHS shell quartet name: SQ_NAI_F_S_C1002_a
   ************************************************************/
  Double I_NAI_F3x_Px_C1002_a = I_NAI_G4x_S_C1002_a+ABX*I_NAI_F3x_S_C1002_a;
  Double I_NAI_F2xy_Px_C1002_a = I_NAI_G3xy_S_C1002_a+ABX*I_NAI_F2xy_S_C1002_a;
  Double I_NAI_F2xz_Px_C1002_a = I_NAI_G3xz_S_C1002_a+ABX*I_NAI_F2xz_S_C1002_a;
  Double I_NAI_Fx2y_Px_C1002_a = I_NAI_G2x2y_S_C1002_a+ABX*I_NAI_Fx2y_S_C1002_a;
  Double I_NAI_Fxyz_Px_C1002_a = I_NAI_G2xyz_S_C1002_a+ABX*I_NAI_Fxyz_S_C1002_a;
  Double I_NAI_Fx2z_Px_C1002_a = I_NAI_G2x2z_S_C1002_a+ABX*I_NAI_Fx2z_S_C1002_a;
  Double I_NAI_F3y_Px_C1002_a = I_NAI_Gx3y_S_C1002_a+ABX*I_NAI_F3y_S_C1002_a;
  Double I_NAI_F2yz_Px_C1002_a = I_NAI_Gx2yz_S_C1002_a+ABX*I_NAI_F2yz_S_C1002_a;
  Double I_NAI_Fy2z_Px_C1002_a = I_NAI_Gxy2z_S_C1002_a+ABX*I_NAI_Fy2z_S_C1002_a;
  Double I_NAI_F3z_Px_C1002_a = I_NAI_Gx3z_S_C1002_a+ABX*I_NAI_F3z_S_C1002_a;
  Double I_NAI_F3x_Py_C1002_a = I_NAI_G3xy_S_C1002_a+ABY*I_NAI_F3x_S_C1002_a;
  Double I_NAI_F2xy_Py_C1002_a = I_NAI_G2x2y_S_C1002_a+ABY*I_NAI_F2xy_S_C1002_a;
  Double I_NAI_F2xz_Py_C1002_a = I_NAI_G2xyz_S_C1002_a+ABY*I_NAI_F2xz_S_C1002_a;
  Double I_NAI_Fx2y_Py_C1002_a = I_NAI_Gx3y_S_C1002_a+ABY*I_NAI_Fx2y_S_C1002_a;
  Double I_NAI_Fxyz_Py_C1002_a = I_NAI_Gx2yz_S_C1002_a+ABY*I_NAI_Fxyz_S_C1002_a;
  Double I_NAI_Fx2z_Py_C1002_a = I_NAI_Gxy2z_S_C1002_a+ABY*I_NAI_Fx2z_S_C1002_a;
  Double I_NAI_F3y_Py_C1002_a = I_NAI_G4y_S_C1002_a+ABY*I_NAI_F3y_S_C1002_a;
  Double I_NAI_F2yz_Py_C1002_a = I_NAI_G3yz_S_C1002_a+ABY*I_NAI_F2yz_S_C1002_a;
  Double I_NAI_Fy2z_Py_C1002_a = I_NAI_G2y2z_S_C1002_a+ABY*I_NAI_Fy2z_S_C1002_a;
  Double I_NAI_F3z_Py_C1002_a = I_NAI_Gy3z_S_C1002_a+ABY*I_NAI_F3z_S_C1002_a;
  Double I_NAI_F3x_Pz_C1002_a = I_NAI_G3xz_S_C1002_a+ABZ*I_NAI_F3x_S_C1002_a;
  Double I_NAI_F2xy_Pz_C1002_a = I_NAI_G2xyz_S_C1002_a+ABZ*I_NAI_F2xy_S_C1002_a;
  Double I_NAI_F2xz_Pz_C1002_a = I_NAI_G2x2z_S_C1002_a+ABZ*I_NAI_F2xz_S_C1002_a;
  Double I_NAI_Fx2y_Pz_C1002_a = I_NAI_Gx2yz_S_C1002_a+ABZ*I_NAI_Fx2y_S_C1002_a;
  Double I_NAI_Fxyz_Pz_C1002_a = I_NAI_Gxy2z_S_C1002_a+ABZ*I_NAI_Fxyz_S_C1002_a;
  Double I_NAI_Fx2z_Pz_C1002_a = I_NAI_Gx3z_S_C1002_a+ABZ*I_NAI_Fx2z_S_C1002_a;
  Double I_NAI_F3y_Pz_C1002_a = I_NAI_G3yz_S_C1002_a+ABZ*I_NAI_F3y_S_C1002_a;
  Double I_NAI_F2yz_Pz_C1002_a = I_NAI_G2y2z_S_C1002_a+ABZ*I_NAI_F2yz_S_C1002_a;
  Double I_NAI_Fy2z_Pz_C1002_a = I_NAI_Gy3z_S_C1002_a+ABZ*I_NAI_Fy2z_S_C1002_a;
  Double I_NAI_F3z_Pz_C1002_a = I_NAI_G4z_S_C1002_a+ABZ*I_NAI_F3z_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C2_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C2_b
   * RHS shell quartet name: SQ_NAI_D_S_C2_b
   ************************************************************/
  Double I_NAI_D2x_Px_C2_b = I_NAI_F3x_S_C2_b+ABX*I_NAI_D2x_S_C2_b;
  Double I_NAI_Dxy_Px_C2_b = I_NAI_F2xy_S_C2_b+ABX*I_NAI_Dxy_S_C2_b;
  Double I_NAI_Dxz_Px_C2_b = I_NAI_F2xz_S_C2_b+ABX*I_NAI_Dxz_S_C2_b;
  Double I_NAI_D2y_Px_C2_b = I_NAI_Fx2y_S_C2_b+ABX*I_NAI_D2y_S_C2_b;
  Double I_NAI_Dyz_Px_C2_b = I_NAI_Fxyz_S_C2_b+ABX*I_NAI_Dyz_S_C2_b;
  Double I_NAI_D2z_Px_C2_b = I_NAI_Fx2z_S_C2_b+ABX*I_NAI_D2z_S_C2_b;
  Double I_NAI_D2x_Py_C2_b = I_NAI_F2xy_S_C2_b+ABY*I_NAI_D2x_S_C2_b;
  Double I_NAI_Dxy_Py_C2_b = I_NAI_Fx2y_S_C2_b+ABY*I_NAI_Dxy_S_C2_b;
  Double I_NAI_Dxz_Py_C2_b = I_NAI_Fxyz_S_C2_b+ABY*I_NAI_Dxz_S_C2_b;
  Double I_NAI_D2y_Py_C2_b = I_NAI_F3y_S_C2_b+ABY*I_NAI_D2y_S_C2_b;
  Double I_NAI_Dyz_Py_C2_b = I_NAI_F2yz_S_C2_b+ABY*I_NAI_Dyz_S_C2_b;
  Double I_NAI_D2z_Py_C2_b = I_NAI_Fy2z_S_C2_b+ABY*I_NAI_D2z_S_C2_b;
  Double I_NAI_D2x_Pz_C2_b = I_NAI_F2xz_S_C2_b+ABZ*I_NAI_D2x_S_C2_b;
  Double I_NAI_Dxy_Pz_C2_b = I_NAI_Fxyz_S_C2_b+ABZ*I_NAI_Dxy_S_C2_b;
  Double I_NAI_Dxz_Pz_C2_b = I_NAI_Fx2z_S_C2_b+ABZ*I_NAI_Dxz_S_C2_b;
  Double I_NAI_D2y_Pz_C2_b = I_NAI_F2yz_S_C2_b+ABZ*I_NAI_D2y_S_C2_b;
  Double I_NAI_Dyz_Pz_C2_b = I_NAI_Fy2z_S_C2_b+ABZ*I_NAI_Dyz_S_C2_b;
  Double I_NAI_D2z_Pz_C2_b = I_NAI_F3z_S_C2_b+ABZ*I_NAI_D2z_S_C2_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C1002_b
   * RHS shell quartet name: SQ_NAI_D_S_C1002_b
   ************************************************************/
  Double I_NAI_D2x_Px_C1002_b = I_NAI_F3x_S_C1002_b+ABX*I_NAI_D2x_S_C1002_b;
  Double I_NAI_Dxy_Px_C1002_b = I_NAI_F2xy_S_C1002_b+ABX*I_NAI_Dxy_S_C1002_b;
  Double I_NAI_Dxz_Px_C1002_b = I_NAI_F2xz_S_C1002_b+ABX*I_NAI_Dxz_S_C1002_b;
  Double I_NAI_D2y_Px_C1002_b = I_NAI_Fx2y_S_C1002_b+ABX*I_NAI_D2y_S_C1002_b;
  Double I_NAI_Dyz_Px_C1002_b = I_NAI_Fxyz_S_C1002_b+ABX*I_NAI_Dyz_S_C1002_b;
  Double I_NAI_D2z_Px_C1002_b = I_NAI_Fx2z_S_C1002_b+ABX*I_NAI_D2z_S_C1002_b;
  Double I_NAI_D2x_Py_C1002_b = I_NAI_F2xy_S_C1002_b+ABY*I_NAI_D2x_S_C1002_b;
  Double I_NAI_Dxy_Py_C1002_b = I_NAI_Fx2y_S_C1002_b+ABY*I_NAI_Dxy_S_C1002_b;
  Double I_NAI_Dxz_Py_C1002_b = I_NAI_Fxyz_S_C1002_b+ABY*I_NAI_Dxz_S_C1002_b;
  Double I_NAI_D2y_Py_C1002_b = I_NAI_F3y_S_C1002_b+ABY*I_NAI_D2y_S_C1002_b;
  Double I_NAI_Dyz_Py_C1002_b = I_NAI_F2yz_S_C1002_b+ABY*I_NAI_Dyz_S_C1002_b;
  Double I_NAI_D2z_Py_C1002_b = I_NAI_Fy2z_S_C1002_b+ABY*I_NAI_D2z_S_C1002_b;
  Double I_NAI_D2x_Pz_C1002_b = I_NAI_F2xz_S_C1002_b+ABZ*I_NAI_D2x_S_C1002_b;
  Double I_NAI_Dxy_Pz_C1002_b = I_NAI_Fxyz_S_C1002_b+ABZ*I_NAI_Dxy_S_C1002_b;
  Double I_NAI_Dxz_Pz_C1002_b = I_NAI_Fx2z_S_C1002_b+ABZ*I_NAI_Dxz_S_C1002_b;
  Double I_NAI_D2y_Pz_C1002_b = I_NAI_F2yz_S_C1002_b+ABZ*I_NAI_D2y_S_C1002_b;
  Double I_NAI_Dyz_Pz_C1002_b = I_NAI_Fy2z_S_C1002_b+ABZ*I_NAI_Dyz_S_C1002_b;
  Double I_NAI_D2z_Pz_C1002_b = I_NAI_F3z_S_C1002_b+ABZ*I_NAI_D2z_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_C1002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 5 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_C1002_b
   * RHS shell quartet name: SQ_NAI_F_S_C1002_b
   ************************************************************/
  Double I_NAI_F3x_Px_C1002_b = I_NAI_G4x_S_C1002_b+ABX*I_NAI_F3x_S_C1002_b;
  Double I_NAI_F2xy_Px_C1002_b = I_NAI_G3xy_S_C1002_b+ABX*I_NAI_F2xy_S_C1002_b;
  Double I_NAI_F2xz_Px_C1002_b = I_NAI_G3xz_S_C1002_b+ABX*I_NAI_F2xz_S_C1002_b;
  Double I_NAI_Fx2y_Px_C1002_b = I_NAI_G2x2y_S_C1002_b+ABX*I_NAI_Fx2y_S_C1002_b;
  Double I_NAI_Fxyz_Px_C1002_b = I_NAI_G2xyz_S_C1002_b+ABX*I_NAI_Fxyz_S_C1002_b;
  Double I_NAI_Fx2z_Px_C1002_b = I_NAI_G2x2z_S_C1002_b+ABX*I_NAI_Fx2z_S_C1002_b;
  Double I_NAI_F3y_Px_C1002_b = I_NAI_Gx3y_S_C1002_b+ABX*I_NAI_F3y_S_C1002_b;
  Double I_NAI_F2yz_Px_C1002_b = I_NAI_Gx2yz_S_C1002_b+ABX*I_NAI_F2yz_S_C1002_b;
  Double I_NAI_Fy2z_Px_C1002_b = I_NAI_Gxy2z_S_C1002_b+ABX*I_NAI_Fy2z_S_C1002_b;
  Double I_NAI_F3z_Px_C1002_b = I_NAI_Gx3z_S_C1002_b+ABX*I_NAI_F3z_S_C1002_b;
  Double I_NAI_F2xy_Py_C1002_b = I_NAI_G2x2y_S_C1002_b+ABY*I_NAI_F2xy_S_C1002_b;
  Double I_NAI_F2xz_Py_C1002_b = I_NAI_G2xyz_S_C1002_b+ABY*I_NAI_F2xz_S_C1002_b;
  Double I_NAI_Fx2y_Py_C1002_b = I_NAI_Gx3y_S_C1002_b+ABY*I_NAI_Fx2y_S_C1002_b;
  Double I_NAI_Fxyz_Py_C1002_b = I_NAI_Gx2yz_S_C1002_b+ABY*I_NAI_Fxyz_S_C1002_b;
  Double I_NAI_Fx2z_Py_C1002_b = I_NAI_Gxy2z_S_C1002_b+ABY*I_NAI_Fx2z_S_C1002_b;
  Double I_NAI_F3y_Py_C1002_b = I_NAI_G4y_S_C1002_b+ABY*I_NAI_F3y_S_C1002_b;
  Double I_NAI_F2yz_Py_C1002_b = I_NAI_G3yz_S_C1002_b+ABY*I_NAI_F2yz_S_C1002_b;
  Double I_NAI_Fy2z_Py_C1002_b = I_NAI_G2y2z_S_C1002_b+ABY*I_NAI_Fy2z_S_C1002_b;
  Double I_NAI_F3z_Py_C1002_b = I_NAI_Gy3z_S_C1002_b+ABY*I_NAI_F3z_S_C1002_b;
  Double I_NAI_F2xz_Pz_C1002_b = I_NAI_G2x2z_S_C1002_b+ABZ*I_NAI_F2xz_S_C1002_b;
  Double I_NAI_Fxyz_Pz_C1002_b = I_NAI_Gxy2z_S_C1002_b+ABZ*I_NAI_Fxyz_S_C1002_b;
  Double I_NAI_Fx2z_Pz_C1002_b = I_NAI_Gx3z_S_C1002_b+ABZ*I_NAI_Fx2z_S_C1002_b;
  Double I_NAI_F2yz_Pz_C1002_b = I_NAI_G2y2z_S_C1002_b+ABZ*I_NAI_F2yz_S_C1002_b;
  Double I_NAI_Fy2z_Pz_C1002_b = I_NAI_Gy3z_S_C1002_b+ABZ*I_NAI_Fy2z_S_C1002_b;
  Double I_NAI_F3z_Pz_C1002_b = I_NAI_G4z_S_C1002_b+ABZ*I_NAI_F3z_S_C1002_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D_C1002_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1002_b
   * RHS shell quartet name: SQ_NAI_D_P_C1002_b
   ************************************************************/
  Double I_NAI_D2x_D2x_C1002_b = I_NAI_F3x_Px_C1002_b+ABX*I_NAI_D2x_Px_C1002_b;
  Double I_NAI_Dxy_D2x_C1002_b = I_NAI_F2xy_Px_C1002_b+ABX*I_NAI_Dxy_Px_C1002_b;
  Double I_NAI_Dxz_D2x_C1002_b = I_NAI_F2xz_Px_C1002_b+ABX*I_NAI_Dxz_Px_C1002_b;
  Double I_NAI_D2y_D2x_C1002_b = I_NAI_Fx2y_Px_C1002_b+ABX*I_NAI_D2y_Px_C1002_b;
  Double I_NAI_Dyz_D2x_C1002_b = I_NAI_Fxyz_Px_C1002_b+ABX*I_NAI_Dyz_Px_C1002_b;
  Double I_NAI_D2z_D2x_C1002_b = I_NAI_Fx2z_Px_C1002_b+ABX*I_NAI_D2z_Px_C1002_b;
  Double I_NAI_D2x_Dxy_C1002_b = I_NAI_F2xy_Px_C1002_b+ABY*I_NAI_D2x_Px_C1002_b;
  Double I_NAI_Dxy_Dxy_C1002_b = I_NAI_Fx2y_Px_C1002_b+ABY*I_NAI_Dxy_Px_C1002_b;
  Double I_NAI_Dxz_Dxy_C1002_b = I_NAI_Fxyz_Px_C1002_b+ABY*I_NAI_Dxz_Px_C1002_b;
  Double I_NAI_D2y_Dxy_C1002_b = I_NAI_F3y_Px_C1002_b+ABY*I_NAI_D2y_Px_C1002_b;
  Double I_NAI_Dyz_Dxy_C1002_b = I_NAI_F2yz_Px_C1002_b+ABY*I_NAI_Dyz_Px_C1002_b;
  Double I_NAI_D2z_Dxy_C1002_b = I_NAI_Fy2z_Px_C1002_b+ABY*I_NAI_D2z_Px_C1002_b;
  Double I_NAI_D2x_Dxz_C1002_b = I_NAI_F2xz_Px_C1002_b+ABZ*I_NAI_D2x_Px_C1002_b;
  Double I_NAI_Dxy_Dxz_C1002_b = I_NAI_Fxyz_Px_C1002_b+ABZ*I_NAI_Dxy_Px_C1002_b;
  Double I_NAI_Dxz_Dxz_C1002_b = I_NAI_Fx2z_Px_C1002_b+ABZ*I_NAI_Dxz_Px_C1002_b;
  Double I_NAI_D2y_Dxz_C1002_b = I_NAI_F2yz_Px_C1002_b+ABZ*I_NAI_D2y_Px_C1002_b;
  Double I_NAI_Dyz_Dxz_C1002_b = I_NAI_Fy2z_Px_C1002_b+ABZ*I_NAI_Dyz_Px_C1002_b;
  Double I_NAI_D2z_Dxz_C1002_b = I_NAI_F3z_Px_C1002_b+ABZ*I_NAI_D2z_Px_C1002_b;
  Double I_NAI_D2x_D2y_C1002_b = I_NAI_F2xy_Py_C1002_b+ABY*I_NAI_D2x_Py_C1002_b;
  Double I_NAI_Dxy_D2y_C1002_b = I_NAI_Fx2y_Py_C1002_b+ABY*I_NAI_Dxy_Py_C1002_b;
  Double I_NAI_Dxz_D2y_C1002_b = I_NAI_Fxyz_Py_C1002_b+ABY*I_NAI_Dxz_Py_C1002_b;
  Double I_NAI_D2y_D2y_C1002_b = I_NAI_F3y_Py_C1002_b+ABY*I_NAI_D2y_Py_C1002_b;
  Double I_NAI_Dyz_D2y_C1002_b = I_NAI_F2yz_Py_C1002_b+ABY*I_NAI_Dyz_Py_C1002_b;
  Double I_NAI_D2z_D2y_C1002_b = I_NAI_Fy2z_Py_C1002_b+ABY*I_NAI_D2z_Py_C1002_b;
  Double I_NAI_D2x_Dyz_C1002_b = I_NAI_F2xz_Py_C1002_b+ABZ*I_NAI_D2x_Py_C1002_b;
  Double I_NAI_Dxy_Dyz_C1002_b = I_NAI_Fxyz_Py_C1002_b+ABZ*I_NAI_Dxy_Py_C1002_b;
  Double I_NAI_Dxz_Dyz_C1002_b = I_NAI_Fx2z_Py_C1002_b+ABZ*I_NAI_Dxz_Py_C1002_b;
  Double I_NAI_D2y_Dyz_C1002_b = I_NAI_F2yz_Py_C1002_b+ABZ*I_NAI_D2y_Py_C1002_b;
  Double I_NAI_Dyz_Dyz_C1002_b = I_NAI_Fy2z_Py_C1002_b+ABZ*I_NAI_Dyz_Py_C1002_b;
  Double I_NAI_D2z_Dyz_C1002_b = I_NAI_F3z_Py_C1002_b+ABZ*I_NAI_D2z_Py_C1002_b;
  Double I_NAI_D2x_D2z_C1002_b = I_NAI_F2xz_Pz_C1002_b+ABZ*I_NAI_D2x_Pz_C1002_b;
  Double I_NAI_Dxy_D2z_C1002_b = I_NAI_Fxyz_Pz_C1002_b+ABZ*I_NAI_Dxy_Pz_C1002_b;
  Double I_NAI_Dxz_D2z_C1002_b = I_NAI_Fx2z_Pz_C1002_b+ABZ*I_NAI_Dxz_Pz_C1002_b;
  Double I_NAI_D2y_D2z_C1002_b = I_NAI_F2yz_Pz_C1002_b+ABZ*I_NAI_D2y_Pz_C1002_b;
  Double I_NAI_Dyz_D2z_C1002_b = I_NAI_Fy2z_Pz_C1002_b+ABZ*I_NAI_Dyz_Pz_C1002_b;
  Double I_NAI_D2z_D2z_C1002_b = I_NAI_F3z_Pz_C1002_b+ABZ*I_NAI_D2z_Pz_C1002_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_S_C2_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C2_a
   * RHS shell quartet name: SQ_NAI_P_S_C2
   ************************************************************/
  abcd[0] = 2.0E0*I_NAI_F3x_S_C2_a-2*I_NAI_Px_S_C2;
  abcd[1] = 2.0E0*I_NAI_F2xy_S_C2_a-1*I_NAI_Py_S_C2;
  abcd[2] = 2.0E0*I_NAI_F2xz_S_C2_a-1*I_NAI_Pz_S_C2;
  abcd[3] = 2.0E0*I_NAI_Fx2y_S_C2_a;
  abcd[4] = 2.0E0*I_NAI_Fxyz_S_C2_a;
  abcd[5] = 2.0E0*I_NAI_Fx2z_S_C2_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1002_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1002_a
   * RHS shell quartet name: SQ_NAI_P_P_C1002
   ************************************************************/
  abcd[6] = 2.0E0*I_NAI_F3x_Px_C1002_a-2*I_NAI_Px_Px_C1002;
  abcd[7] = 2.0E0*I_NAI_F2xy_Px_C1002_a-1*I_NAI_Py_Px_C1002;
  abcd[8] = 2.0E0*I_NAI_F2xz_Px_C1002_a-1*I_NAI_Pz_Px_C1002;
  abcd[9] = 2.0E0*I_NAI_Fx2y_Px_C1002_a;
  abcd[10] = 2.0E0*I_NAI_Fxyz_Px_C1002_a;
  abcd[11] = 2.0E0*I_NAI_Fx2z_Px_C1002_a;
  abcd[12] = 2.0E0*I_NAI_F3x_Py_C1002_a-2*I_NAI_Px_Py_C1002;
  abcd[13] = 2.0E0*I_NAI_F2xy_Py_C1002_a-1*I_NAI_Py_Py_C1002;
  abcd[14] = 2.0E0*I_NAI_F2xz_Py_C1002_a-1*I_NAI_Pz_Py_C1002;
  abcd[15] = 2.0E0*I_NAI_Fx2y_Py_C1002_a;
  abcd[16] = 2.0E0*I_NAI_Fxyz_Py_C1002_a;
  abcd[17] = 2.0E0*I_NAI_Fx2z_Py_C1002_a;
  abcd[18] = 2.0E0*I_NAI_F3x_Pz_C1002_a-2*I_NAI_Px_Pz_C1002;
  abcd[19] = 2.0E0*I_NAI_F2xy_Pz_C1002_a-1*I_NAI_Py_Pz_C1002;
  abcd[20] = 2.0E0*I_NAI_F2xz_Pz_C1002_a-1*I_NAI_Pz_Pz_C1002;
  abcd[21] = 2.0E0*I_NAI_Fx2y_Pz_C1002_a;
  abcd[22] = 2.0E0*I_NAI_Fxyz_Pz_C1002_a;
  abcd[23] = 2.0E0*I_NAI_Fx2z_Pz_C1002_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_S_C2_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C2_a
   * RHS shell quartet name: SQ_NAI_P_S_C2
   ************************************************************/
  abcd[24] = 2.0E0*I_NAI_F2xy_S_C2_a;
  abcd[25] = 2.0E0*I_NAI_Fx2y_S_C2_a-1*I_NAI_Px_S_C2;
  abcd[26] = 2.0E0*I_NAI_Fxyz_S_C2_a;
  abcd[27] = 2.0E0*I_NAI_F3y_S_C2_a-2*I_NAI_Py_S_C2;
  abcd[28] = 2.0E0*I_NAI_F2yz_S_C2_a-1*I_NAI_Pz_S_C2;
  abcd[29] = 2.0E0*I_NAI_Fy2z_S_C2_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1002_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1002_a
   * RHS shell quartet name: SQ_NAI_P_P_C1002
   ************************************************************/
  abcd[30] = 2.0E0*I_NAI_F2xy_Px_C1002_a;
  abcd[31] = 2.0E0*I_NAI_Fx2y_Px_C1002_a-1*I_NAI_Px_Px_C1002;
  abcd[32] = 2.0E0*I_NAI_Fxyz_Px_C1002_a;
  abcd[33] = 2.0E0*I_NAI_F3y_Px_C1002_a-2*I_NAI_Py_Px_C1002;
  abcd[34] = 2.0E0*I_NAI_F2yz_Px_C1002_a-1*I_NAI_Pz_Px_C1002;
  abcd[35] = 2.0E0*I_NAI_Fy2z_Px_C1002_a;
  abcd[36] = 2.0E0*I_NAI_F2xy_Py_C1002_a;
  abcd[37] = 2.0E0*I_NAI_Fx2y_Py_C1002_a-1*I_NAI_Px_Py_C1002;
  abcd[38] = 2.0E0*I_NAI_Fxyz_Py_C1002_a;
  abcd[39] = 2.0E0*I_NAI_F3y_Py_C1002_a-2*I_NAI_Py_Py_C1002;
  abcd[40] = 2.0E0*I_NAI_F2yz_Py_C1002_a-1*I_NAI_Pz_Py_C1002;
  abcd[41] = 2.0E0*I_NAI_Fy2z_Py_C1002_a;
  abcd[42] = 2.0E0*I_NAI_F2xy_Pz_C1002_a;
  abcd[43] = 2.0E0*I_NAI_Fx2y_Pz_C1002_a-1*I_NAI_Px_Pz_C1002;
  abcd[44] = 2.0E0*I_NAI_Fxyz_Pz_C1002_a;
  abcd[45] = 2.0E0*I_NAI_F3y_Pz_C1002_a-2*I_NAI_Py_Pz_C1002;
  abcd[46] = 2.0E0*I_NAI_F2yz_Pz_C1002_a-1*I_NAI_Pz_Pz_C1002;
  abcd[47] = 2.0E0*I_NAI_Fy2z_Pz_C1002_a;

  /************************************************************
   * shell quartet name: SQ_NAI_D_S_C2_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S_C2_a
   * RHS shell quartet name: SQ_NAI_P_S_C2
   ************************************************************/
  abcd[48] = 2.0E0*I_NAI_F2xz_S_C2_a;
  abcd[49] = 2.0E0*I_NAI_Fxyz_S_C2_a;
  abcd[50] = 2.0E0*I_NAI_Fx2z_S_C2_a-1*I_NAI_Px_S_C2;
  abcd[51] = 2.0E0*I_NAI_F2yz_S_C2_a;
  abcd[52] = 2.0E0*I_NAI_Fy2z_S_C2_a-1*I_NAI_Py_S_C2;
  abcd[53] = 2.0E0*I_NAI_F3z_S_C2_a-2*I_NAI_Pz_S_C2;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1002_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P_C1002_a
   * RHS shell quartet name: SQ_NAI_P_P_C1002
   ************************************************************/
  abcd[54] = 2.0E0*I_NAI_F2xz_Px_C1002_a;
  abcd[55] = 2.0E0*I_NAI_Fxyz_Px_C1002_a;
  abcd[56] = 2.0E0*I_NAI_Fx2z_Px_C1002_a-1*I_NAI_Px_Px_C1002;
  abcd[57] = 2.0E0*I_NAI_F2yz_Px_C1002_a;
  abcd[58] = 2.0E0*I_NAI_Fy2z_Px_C1002_a-1*I_NAI_Py_Px_C1002;
  abcd[59] = 2.0E0*I_NAI_F3z_Px_C1002_a-2*I_NAI_Pz_Px_C1002;
  abcd[60] = 2.0E0*I_NAI_F2xz_Py_C1002_a;
  abcd[61] = 2.0E0*I_NAI_Fxyz_Py_C1002_a;
  abcd[62] = 2.0E0*I_NAI_Fx2z_Py_C1002_a-1*I_NAI_Px_Py_C1002;
  abcd[63] = 2.0E0*I_NAI_F2yz_Py_C1002_a;
  abcd[64] = 2.0E0*I_NAI_Fy2z_Py_C1002_a-1*I_NAI_Py_Py_C1002;
  abcd[65] = 2.0E0*I_NAI_F3z_Py_C1002_a-2*I_NAI_Pz_Py_C1002;
  abcd[66] = 2.0E0*I_NAI_F2xz_Pz_C1002_a;
  abcd[67] = 2.0E0*I_NAI_Fxyz_Pz_C1002_a;
  abcd[68] = 2.0E0*I_NAI_Fx2z_Pz_C1002_a-1*I_NAI_Px_Pz_C1002;
  abcd[69] = 2.0E0*I_NAI_F2yz_Pz_C1002_a;
  abcd[70] = 2.0E0*I_NAI_Fy2z_Pz_C1002_a-1*I_NAI_Py_Pz_C1002;
  abcd[71] = 2.0E0*I_NAI_F3z_Pz_C1002_a-2*I_NAI_Pz_Pz_C1002;

  /************************************************************
   * shell quartet name: SQ_NAI_D_S_C2_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C2_b
   ************************************************************/
  abcd[72] = 2.0E0*I_NAI_D2x_Px_C2_b;
  abcd[73] = 2.0E0*I_NAI_Dxy_Px_C2_b;
  abcd[74] = 2.0E0*I_NAI_Dxz_Px_C2_b;
  abcd[75] = 2.0E0*I_NAI_D2y_Px_C2_b;
  abcd[76] = 2.0E0*I_NAI_Dyz_Px_C2_b;
  abcd[77] = 2.0E0*I_NAI_D2z_Px_C2_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1002_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1002_b
   * RHS shell quartet name: SQ_NAI_D_S_C1002
   ************************************************************/
  abcd[78] = 2.0E0*I_NAI_D2x_D2x_C1002_b-1*I_NAI_D2x_S_C1002;
  abcd[79] = 2.0E0*I_NAI_Dxy_D2x_C1002_b-1*I_NAI_Dxy_S_C1002;
  abcd[80] = 2.0E0*I_NAI_Dxz_D2x_C1002_b-1*I_NAI_Dxz_S_C1002;
  abcd[81] = 2.0E0*I_NAI_D2y_D2x_C1002_b-1*I_NAI_D2y_S_C1002;
  abcd[82] = 2.0E0*I_NAI_Dyz_D2x_C1002_b-1*I_NAI_Dyz_S_C1002;
  abcd[83] = 2.0E0*I_NAI_D2z_D2x_C1002_b-1*I_NAI_D2z_S_C1002;
  abcd[84] = 2.0E0*I_NAI_D2x_Dxy_C1002_b;
  abcd[85] = 2.0E0*I_NAI_Dxy_Dxy_C1002_b;
  abcd[86] = 2.0E0*I_NAI_Dxz_Dxy_C1002_b;
  abcd[87] = 2.0E0*I_NAI_D2y_Dxy_C1002_b;
  abcd[88] = 2.0E0*I_NAI_Dyz_Dxy_C1002_b;
  abcd[89] = 2.0E0*I_NAI_D2z_Dxy_C1002_b;
  abcd[90] = 2.0E0*I_NAI_D2x_Dxz_C1002_b;
  abcd[91] = 2.0E0*I_NAI_Dxy_Dxz_C1002_b;
  abcd[92] = 2.0E0*I_NAI_Dxz_Dxz_C1002_b;
  abcd[93] = 2.0E0*I_NAI_D2y_Dxz_C1002_b;
  abcd[94] = 2.0E0*I_NAI_Dyz_Dxz_C1002_b;
  abcd[95] = 2.0E0*I_NAI_D2z_Dxz_C1002_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_S_C2_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C2_b
   ************************************************************/
  abcd[96] = 2.0E0*I_NAI_D2x_Py_C2_b;
  abcd[97] = 2.0E0*I_NAI_Dxy_Py_C2_b;
  abcd[98] = 2.0E0*I_NAI_Dxz_Py_C2_b;
  abcd[99] = 2.0E0*I_NAI_D2y_Py_C2_b;
  abcd[100] = 2.0E0*I_NAI_Dyz_Py_C2_b;
  abcd[101] = 2.0E0*I_NAI_D2z_Py_C2_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1002_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1002_b
   * RHS shell quartet name: SQ_NAI_D_S_C1002
   ************************************************************/
  abcd[102] = 2.0E0*I_NAI_D2x_Dxy_C1002_b;
  abcd[103] = 2.0E0*I_NAI_Dxy_Dxy_C1002_b;
  abcd[104] = 2.0E0*I_NAI_Dxz_Dxy_C1002_b;
  abcd[105] = 2.0E0*I_NAI_D2y_Dxy_C1002_b;
  abcd[106] = 2.0E0*I_NAI_Dyz_Dxy_C1002_b;
  abcd[107] = 2.0E0*I_NAI_D2z_Dxy_C1002_b;
  abcd[108] = 2.0E0*I_NAI_D2x_D2y_C1002_b-1*I_NAI_D2x_S_C1002;
  abcd[109] = 2.0E0*I_NAI_Dxy_D2y_C1002_b-1*I_NAI_Dxy_S_C1002;
  abcd[110] = 2.0E0*I_NAI_Dxz_D2y_C1002_b-1*I_NAI_Dxz_S_C1002;
  abcd[111] = 2.0E0*I_NAI_D2y_D2y_C1002_b-1*I_NAI_D2y_S_C1002;
  abcd[112] = 2.0E0*I_NAI_Dyz_D2y_C1002_b-1*I_NAI_Dyz_S_C1002;
  abcd[113] = 2.0E0*I_NAI_D2z_D2y_C1002_b-1*I_NAI_D2z_S_C1002;
  abcd[114] = 2.0E0*I_NAI_D2x_Dyz_C1002_b;
  abcd[115] = 2.0E0*I_NAI_Dxy_Dyz_C1002_b;
  abcd[116] = 2.0E0*I_NAI_Dxz_Dyz_C1002_b;
  abcd[117] = 2.0E0*I_NAI_D2y_Dyz_C1002_b;
  abcd[118] = 2.0E0*I_NAI_Dyz_Dyz_C1002_b;
  abcd[119] = 2.0E0*I_NAI_D2z_Dyz_C1002_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_S_C2_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_P_C2_b
   ************************************************************/
  abcd[120] = 2.0E0*I_NAI_D2x_Pz_C2_b;
  abcd[121] = 2.0E0*I_NAI_Dxy_Pz_C2_b;
  abcd[122] = 2.0E0*I_NAI_Dxz_Pz_C2_b;
  abcd[123] = 2.0E0*I_NAI_D2y_Pz_C2_b;
  abcd[124] = 2.0E0*I_NAI_Dyz_Pz_C2_b;
  abcd[125] = 2.0E0*I_NAI_D2z_Pz_C2_b;

  /************************************************************
   * shell quartet name: SQ_NAI_D_P_C1002_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_D_D_C1002_b
   * RHS shell quartet name: SQ_NAI_D_S_C1002
   ************************************************************/
  abcd[126] = 2.0E0*I_NAI_D2x_Dxz_C1002_b;
  abcd[127] = 2.0E0*I_NAI_Dxy_Dxz_C1002_b;
  abcd[128] = 2.0E0*I_NAI_Dxz_Dxz_C1002_b;
  abcd[129] = 2.0E0*I_NAI_D2y_Dxz_C1002_b;
  abcd[130] = 2.0E0*I_NAI_Dyz_Dxz_C1002_b;
  abcd[131] = 2.0E0*I_NAI_D2z_Dxz_C1002_b;
  abcd[132] = 2.0E0*I_NAI_D2x_Dyz_C1002_b;
  abcd[133] = 2.0E0*I_NAI_Dxy_Dyz_C1002_b;
  abcd[134] = 2.0E0*I_NAI_Dxz_Dyz_C1002_b;
  abcd[135] = 2.0E0*I_NAI_D2y_Dyz_C1002_b;
  abcd[136] = 2.0E0*I_NAI_Dyz_Dyz_C1002_b;
  abcd[137] = 2.0E0*I_NAI_D2z_Dyz_C1002_b;
  abcd[138] = 2.0E0*I_NAI_D2x_D2z_C1002_b-1*I_NAI_D2x_S_C1002;
  abcd[139] = 2.0E0*I_NAI_Dxy_D2z_C1002_b-1*I_NAI_Dxy_S_C1002;
  abcd[140] = 2.0E0*I_NAI_Dxz_D2z_C1002_b-1*I_NAI_Dxz_S_C1002;
  abcd[141] = 2.0E0*I_NAI_D2y_D2z_C1002_b-1*I_NAI_D2y_S_C1002;
  abcd[142] = 2.0E0*I_NAI_Dyz_D2z_C1002_b-1*I_NAI_Dyz_S_C1002;
  abcd[143] = 2.0E0*I_NAI_D2z_D2z_C1002_b-1*I_NAI_D2z_S_C1002;
}
