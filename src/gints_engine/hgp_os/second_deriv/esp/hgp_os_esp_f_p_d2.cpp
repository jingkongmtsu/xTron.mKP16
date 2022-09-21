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
// BRA1 as redundant position, total RHS integrals evaluated as: 2492
// BRA2 as redundant position, total RHS integrals evaluated as: 2136
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

void hgp_os_esp_f_p_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_I6x_S_aa = 0.0E0;
    Double I_ESP_I5xy_S_aa = 0.0E0;
    Double I_ESP_I5xz_S_aa = 0.0E0;
    Double I_ESP_I4x2y_S_aa = 0.0E0;
    Double I_ESP_I4xyz_S_aa = 0.0E0;
    Double I_ESP_I4x2z_S_aa = 0.0E0;
    Double I_ESP_I3x3y_S_aa = 0.0E0;
    Double I_ESP_I3x2yz_S_aa = 0.0E0;
    Double I_ESP_I3xy2z_S_aa = 0.0E0;
    Double I_ESP_I3x3z_S_aa = 0.0E0;
    Double I_ESP_I2x4y_S_aa = 0.0E0;
    Double I_ESP_I2x3yz_S_aa = 0.0E0;
    Double I_ESP_I2x2y2z_S_aa = 0.0E0;
    Double I_ESP_I2xy3z_S_aa = 0.0E0;
    Double I_ESP_I2x4z_S_aa = 0.0E0;
    Double I_ESP_Ix5y_S_aa = 0.0E0;
    Double I_ESP_Ix4yz_S_aa = 0.0E0;
    Double I_ESP_Ix3y2z_S_aa = 0.0E0;
    Double I_ESP_Ix2y3z_S_aa = 0.0E0;
    Double I_ESP_Ixy4z_S_aa = 0.0E0;
    Double I_ESP_Ix5z_S_aa = 0.0E0;
    Double I_ESP_I6y_S_aa = 0.0E0;
    Double I_ESP_I5yz_S_aa = 0.0E0;
    Double I_ESP_I4y2z_S_aa = 0.0E0;
    Double I_ESP_I3y3z_S_aa = 0.0E0;
    Double I_ESP_I2y4z_S_aa = 0.0E0;
    Double I_ESP_Iy5z_S_aa = 0.0E0;
    Double I_ESP_I6z_S_aa = 0.0E0;
    Double I_ESP_H5x_S_aa = 0.0E0;
    Double I_ESP_H4xy_S_aa = 0.0E0;
    Double I_ESP_H4xz_S_aa = 0.0E0;
    Double I_ESP_H3x2y_S_aa = 0.0E0;
    Double I_ESP_H3xyz_S_aa = 0.0E0;
    Double I_ESP_H3x2z_S_aa = 0.0E0;
    Double I_ESP_H2x3y_S_aa = 0.0E0;
    Double I_ESP_H2x2yz_S_aa = 0.0E0;
    Double I_ESP_H2xy2z_S_aa = 0.0E0;
    Double I_ESP_H2x3z_S_aa = 0.0E0;
    Double I_ESP_Hx4y_S_aa = 0.0E0;
    Double I_ESP_Hx3yz_S_aa = 0.0E0;
    Double I_ESP_Hx2y2z_S_aa = 0.0E0;
    Double I_ESP_Hxy3z_S_aa = 0.0E0;
    Double I_ESP_Hx4z_S_aa = 0.0E0;
    Double I_ESP_H5y_S_aa = 0.0E0;
    Double I_ESP_H4yz_S_aa = 0.0E0;
    Double I_ESP_H3y2z_S_aa = 0.0E0;
    Double I_ESP_H2y3z_S_aa = 0.0E0;
    Double I_ESP_Hy4z_S_aa = 0.0E0;
    Double I_ESP_H5z_S_aa = 0.0E0;
    Double I_ESP_G4x_S_a = 0.0E0;
    Double I_ESP_G3xy_S_a = 0.0E0;
    Double I_ESP_G3xz_S_a = 0.0E0;
    Double I_ESP_G2x2y_S_a = 0.0E0;
    Double I_ESP_G2xyz_S_a = 0.0E0;
    Double I_ESP_G2x2z_S_a = 0.0E0;
    Double I_ESP_Gx3y_S_a = 0.0E0;
    Double I_ESP_Gx2yz_S_a = 0.0E0;
    Double I_ESP_Gxy2z_S_a = 0.0E0;
    Double I_ESP_Gx3z_S_a = 0.0E0;
    Double I_ESP_G4y_S_a = 0.0E0;
    Double I_ESP_G3yz_S_a = 0.0E0;
    Double I_ESP_G2y2z_S_a = 0.0E0;
    Double I_ESP_Gy3z_S_a = 0.0E0;
    Double I_ESP_G4z_S_a = 0.0E0;
    Double I_ESP_F3x_S_a = 0.0E0;
    Double I_ESP_F2xy_S_a = 0.0E0;
    Double I_ESP_F2xz_S_a = 0.0E0;
    Double I_ESP_Fx2y_S_a = 0.0E0;
    Double I_ESP_Fxyz_S_a = 0.0E0;
    Double I_ESP_Fx2z_S_a = 0.0E0;
    Double I_ESP_F3y_S_a = 0.0E0;
    Double I_ESP_F2yz_S_a = 0.0E0;
    Double I_ESP_Fy2z_S_a = 0.0E0;
    Double I_ESP_F3z_S_a = 0.0E0;
    Double I_ESP_D2x_S = 0.0E0;
    Double I_ESP_Dxy_S = 0.0E0;
    Double I_ESP_Dxz_S = 0.0E0;
    Double I_ESP_D2y_S = 0.0E0;
    Double I_ESP_Dyz_S = 0.0E0;
    Double I_ESP_D2z_S = 0.0E0;
    Double I_ESP_Px_S = 0.0E0;
    Double I_ESP_Py_S = 0.0E0;
    Double I_ESP_Pz_S = 0.0E0;

    for(UInt ip2=0; ip2<inp2; ip2++) {
      Double ic2   = icoe[ip2];
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
      Double PRX   = PX - R[iGrid*3  ];
      Double PRY   = PY - R[iGrid*3+1];
      Double PRZ   = PZ - R[iGrid*3+2];
      Double PR2   = PRX*PRX+PRY*PRY+PRZ*PRZ;
      Double u     = rho*PR2;
      Double squ   = sqrt(u);
      Double prefactor = ic2*fbra;

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

      Double I_ESP_S_S_vrr  = 0.0E0;
      Double I_ESP_S_S_M1_vrr  = 0.0E0;
      Double I_ESP_S_S_M2_vrr  = 0.0E0;
      Double I_ESP_S_S_M3_vrr  = 0.0E0;
      Double I_ESP_S_S_M4_vrr  = 0.0E0;
      Double I_ESP_S_S_M5_vrr  = 0.0E0;
      Double I_ESP_S_S_M6_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ESP_S_S_vrr_d  = 0.0E0;
      double I_ESP_S_S_M1_vrr_d  = 0.0E0;
      double I_ESP_S_S_M2_vrr_d  = 0.0E0;
      double I_ESP_S_S_M3_vrr_d  = 0.0E0;
      double I_ESP_S_S_M4_vrr_d  = 0.0E0;
      double I_ESP_S_S_M5_vrr_d  = 0.0E0;
      double I_ESP_S_S_M6_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER47;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER45*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER43*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER41*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER39*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER37*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER35*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER33*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER31*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER29*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER27*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER25*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER23*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER21*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER19*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER17*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = 1.0E0+u2*ONEOVER15*I_ESP_S_S_M6_vrr;
        I_ESP_S_S_M6_vrr = ONEOVER13*I_ESP_S_S_M6_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M6_vrr  = f*I_ESP_S_S_M6_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ESP_S_S_M5_vrr  = ONEOVER11*(u2*I_ESP_S_S_M6_vrr+f);
        I_ESP_S_S_M4_vrr  = ONEOVER9*(u2*I_ESP_S_S_M5_vrr+f);
        I_ESP_S_S_M3_vrr  = ONEOVER7*(u2*I_ESP_S_S_M4_vrr+f);
        I_ESP_S_S_M2_vrr  = ONEOVER5*(u2*I_ESP_S_S_M3_vrr+f);
        I_ESP_S_S_M1_vrr  = ONEOVER3*(u2*I_ESP_S_S_M2_vrr+f);
        I_ESP_S_S_vrr  = ONEOVER1*(u2*I_ESP_S_S_M1_vrr+f);

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
          I_ESP_S_S_vrr_d = fac_d*sqrho_d*TWOOVERSQRTPI;
        }else{
          I_ESP_S_S_vrr_d = (fac_d*sqrho_d/squ_d)*erf(squ_d);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        double oneO2u  = 0.5E0/u_d;
        double eu      = exp(-u_d);
        double f       = TWOOVERSQRTPI*fac_d*sqrho_d*eu;
        I_ESP_S_S_M1_vrr_d = oneO2u*(1.0E0*I_ESP_S_S_vrr_d-f);
        I_ESP_S_S_M2_vrr_d = oneO2u*(3.0E0*I_ESP_S_S_M1_vrr_d-f);
        I_ESP_S_S_M3_vrr_d = oneO2u*(5.0E0*I_ESP_S_S_M2_vrr_d-f);
        I_ESP_S_S_M4_vrr_d = oneO2u*(7.0E0*I_ESP_S_S_M3_vrr_d-f);
        I_ESP_S_S_M5_vrr_d = oneO2u*(9.0E0*I_ESP_S_S_M4_vrr_d-f);
        I_ESP_S_S_M6_vrr_d = oneO2u*(11.0E0*I_ESP_S_S_M5_vrr_d-f);

        // write the double result back to the float var
        I_ESP_S_S_vrr = static_cast<Double>(I_ESP_S_S_vrr_d);
        I_ESP_S_S_M1_vrr = static_cast<Double>(I_ESP_S_S_M1_vrr_d);
        I_ESP_S_S_M2_vrr = static_cast<Double>(I_ESP_S_S_M2_vrr_d);
        I_ESP_S_S_M3_vrr = static_cast<Double>(I_ESP_S_S_M3_vrr_d);
        I_ESP_S_S_M4_vrr = static_cast<Double>(I_ESP_S_S_M4_vrr_d);
        I_ESP_S_S_M5_vrr = static_cast<Double>(I_ESP_S_S_M5_vrr_d);
        I_ESP_S_S_M6_vrr = static_cast<Double>(I_ESP_S_S_M6_vrr_d);

#else

        // use erf function to get (SS|SS)^{0}
        if (fabs(u)<THRESHOLD_MATH) {
          I_ESP_S_S_vrr = prefactor*sqrho*TWOOVERSQRTPI;
        }else{
          I_ESP_S_S_vrr = (prefactor*sqrho/squ)*erf(squ);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        Double oneO2u = 0.5E0/u;
        Double eu     = exp(-u);
        Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M1_vrr = oneO2u*(1.0E0*I_ESP_S_S_vrr-f);
        I_ESP_S_S_M2_vrr = oneO2u*(3.0E0*I_ESP_S_S_M1_vrr-f);
        I_ESP_S_S_M3_vrr = oneO2u*(5.0E0*I_ESP_S_S_M2_vrr-f);
        I_ESP_S_S_M4_vrr = oneO2u*(7.0E0*I_ESP_S_S_M3_vrr-f);
        I_ESP_S_S_M5_vrr = oneO2u*(9.0E0*I_ESP_S_S_M4_vrr-f);
        I_ESP_S_S_M6_vrr = oneO2u*(11.0E0*I_ESP_S_S_M5_vrr-f);

#endif

      }


        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M5
         * RHS shell quartet name: SQ_ESP_S_S_M6
         ************************************************************/
        Double I_ESP_Px_S_M5_vrr = PAX*I_ESP_S_S_M5_vrr-PRX*I_ESP_S_S_M6_vrr;
        Double I_ESP_Py_S_M5_vrr = PAY*I_ESP_S_S_M5_vrr-PRY*I_ESP_S_S_M6_vrr;
        Double I_ESP_Pz_S_M5_vrr = PAZ*I_ESP_S_S_M5_vrr-PRZ*I_ESP_S_S_M6_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M5
         ************************************************************/
        Double I_ESP_Px_S_M4_vrr = PAX*I_ESP_S_S_M4_vrr-PRX*I_ESP_S_S_M5_vrr;
        Double I_ESP_Py_S_M4_vrr = PAY*I_ESP_S_S_M4_vrr-PRY*I_ESP_S_S_M5_vrr;
        Double I_ESP_Pz_S_M4_vrr = PAZ*I_ESP_S_S_M4_vrr-PRZ*I_ESP_S_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M4
         * RHS shell quartet name: SQ_ESP_P_S_M5
         * RHS shell quartet name: SQ_ESP_S_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M5
         ************************************************************/
        Double I_ESP_D2x_S_M4_vrr = PAX*I_ESP_Px_S_M4_vrr-PRX*I_ESP_Px_S_M5_vrr+oned2z*I_ESP_S_S_M4_vrr-oned2z*I_ESP_S_S_M5_vrr;
        Double I_ESP_D2y_S_M4_vrr = PAY*I_ESP_Py_S_M4_vrr-PRY*I_ESP_Py_S_M5_vrr+oned2z*I_ESP_S_S_M4_vrr-oned2z*I_ESP_S_S_M5_vrr;
        Double I_ESP_D2z_S_M4_vrr = PAZ*I_ESP_Pz_S_M4_vrr-PRZ*I_ESP_Pz_S_M5_vrr+oned2z*I_ESP_S_S_M4_vrr-oned2z*I_ESP_S_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M4
         ************************************************************/
        Double I_ESP_Px_S_M3_vrr = PAX*I_ESP_S_S_M3_vrr-PRX*I_ESP_S_S_M4_vrr;
        Double I_ESP_Py_S_M3_vrr = PAY*I_ESP_S_S_M3_vrr-PRY*I_ESP_S_S_M4_vrr;
        Double I_ESP_Pz_S_M3_vrr = PAZ*I_ESP_S_S_M3_vrr-PRZ*I_ESP_S_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M4
         ************************************************************/
        Double I_ESP_D2x_S_M3_vrr = PAX*I_ESP_Px_S_M3_vrr-PRX*I_ESP_Px_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;
        Double I_ESP_D2y_S_M3_vrr = PAY*I_ESP_Py_S_M3_vrr-PRY*I_ESP_Py_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;
        Double I_ESP_D2z_S_M3_vrr = PAZ*I_ESP_Pz_S_M3_vrr-PRZ*I_ESP_Pz_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 6 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M3
         * RHS shell quartet name: SQ_ESP_D_S_M4
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M4
         ************************************************************/
        Double I_ESP_F3x_S_M3_vrr = PAX*I_ESP_D2x_S_M3_vrr-PRX*I_ESP_D2x_S_M4_vrr+2*oned2z*I_ESP_Px_S_M3_vrr-2*oned2z*I_ESP_Px_S_M4_vrr;
        Double I_ESP_F2xy_S_M3_vrr = PAY*I_ESP_D2x_S_M3_vrr-PRY*I_ESP_D2x_S_M4_vrr;
        Double I_ESP_F3y_S_M3_vrr = PAY*I_ESP_D2y_S_M3_vrr-PRY*I_ESP_D2y_S_M4_vrr+2*oned2z*I_ESP_Py_S_M3_vrr-2*oned2z*I_ESP_Py_S_M4_vrr;
        Double I_ESP_F3z_S_M3_vrr = PAZ*I_ESP_D2z_S_M3_vrr-PRZ*I_ESP_D2z_S_M4_vrr+2*oned2z*I_ESP_Pz_S_M3_vrr-2*oned2z*I_ESP_Pz_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M3
         ************************************************************/
        Double I_ESP_Px_S_M2_vrr = PAX*I_ESP_S_S_M2_vrr-PRX*I_ESP_S_S_M3_vrr;
        Double I_ESP_Py_S_M2_vrr = PAY*I_ESP_S_S_M2_vrr-PRY*I_ESP_S_S_M3_vrr;
        Double I_ESP_Pz_S_M2_vrr = PAZ*I_ESP_S_S_M2_vrr-PRZ*I_ESP_S_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M3
         ************************************************************/
        Double I_ESP_D2x_S_M2_vrr = PAX*I_ESP_Px_S_M2_vrr-PRX*I_ESP_Px_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;
        Double I_ESP_D2y_S_M2_vrr = PAY*I_ESP_Py_S_M2_vrr-PRY*I_ESP_Py_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;
        Double I_ESP_D2z_S_M2_vrr = PAZ*I_ESP_Pz_S_M2_vrr-PRZ*I_ESP_Pz_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 4 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M3
         ************************************************************/
        Double I_ESP_F3x_S_M2_vrr = PAX*I_ESP_D2x_S_M2_vrr-PRX*I_ESP_D2x_S_M3_vrr+2*oned2z*I_ESP_Px_S_M2_vrr-2*oned2z*I_ESP_Px_S_M3_vrr;
        Double I_ESP_F2xy_S_M2_vrr = PAY*I_ESP_D2x_S_M2_vrr-PRY*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_F2xz_S_M2_vrr = PAZ*I_ESP_D2x_S_M2_vrr-PRZ*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_F3y_S_M2_vrr = PAY*I_ESP_D2y_S_M2_vrr-PRY*I_ESP_D2y_S_M3_vrr+2*oned2z*I_ESP_Py_S_M2_vrr-2*oned2z*I_ESP_Py_S_M3_vrr;
        Double I_ESP_F2yz_S_M2_vrr = PAZ*I_ESP_D2y_S_M2_vrr-PRZ*I_ESP_D2y_S_M3_vrr;
        Double I_ESP_F3z_S_M2_vrr = PAZ*I_ESP_D2z_S_M2_vrr-PRZ*I_ESP_D2z_S_M3_vrr+2*oned2z*I_ESP_Pz_S_M2_vrr-2*oned2z*I_ESP_Pz_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 5 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M2
         * RHS shell quartet name: SQ_ESP_F_S_M3
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M3
         ************************************************************/
        Double I_ESP_G4x_S_M2_vrr = PAX*I_ESP_F3x_S_M2_vrr-PRX*I_ESP_F3x_S_M3_vrr+3*oned2z*I_ESP_D2x_S_M2_vrr-3*oned2z*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_G3xy_S_M2_vrr = PAY*I_ESP_F3x_S_M2_vrr-PRY*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_G3xz_S_M2_vrr = PAZ*I_ESP_F3x_S_M2_vrr-PRZ*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_G2x2y_S_M2_vrr = PAY*I_ESP_F2xy_S_M2_vrr-PRY*I_ESP_F2xy_S_M3_vrr+oned2z*I_ESP_D2x_S_M2_vrr-oned2z*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_Gx3y_S_M2_vrr = PAX*I_ESP_F3y_S_M2_vrr-PRX*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_Gx3z_S_M2_vrr = PAX*I_ESP_F3z_S_M2_vrr-PRX*I_ESP_F3z_S_M3_vrr;
        Double I_ESP_G4y_S_M2_vrr = PAY*I_ESP_F3y_S_M2_vrr-PRY*I_ESP_F3y_S_M3_vrr+3*oned2z*I_ESP_D2y_S_M2_vrr-3*oned2z*I_ESP_D2y_S_M3_vrr;
        Double I_ESP_G3yz_S_M2_vrr = PAZ*I_ESP_F3y_S_M2_vrr-PRZ*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_Gy3z_S_M2_vrr = PAY*I_ESP_F3z_S_M2_vrr-PRY*I_ESP_F3z_S_M3_vrr;
        Double I_ESP_G4z_S_M2_vrr = PAZ*I_ESP_F3z_S_M2_vrr-PRZ*I_ESP_F3z_S_M3_vrr+3*oned2z*I_ESP_D2z_S_M2_vrr-3*oned2z*I_ESP_D2z_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_Px_S_M1_vrr = PAX*I_ESP_S_S_M1_vrr-PRX*I_ESP_S_S_M2_vrr;
        Double I_ESP_Py_S_M1_vrr = PAY*I_ESP_S_S_M1_vrr-PRY*I_ESP_S_S_M2_vrr;
        Double I_ESP_Pz_S_M1_vrr = PAZ*I_ESP_S_S_M1_vrr-PRZ*I_ESP_S_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_D2x_S_M1_vrr = PAX*I_ESP_Px_S_M1_vrr-PRX*I_ESP_Px_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
        Double I_ESP_Dxy_S_M1_vrr = PAY*I_ESP_Px_S_M1_vrr-PRY*I_ESP_Px_S_M2_vrr;
        Double I_ESP_D2y_S_M1_vrr = PAY*I_ESP_Py_S_M1_vrr-PRY*I_ESP_Py_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
        Double I_ESP_D2z_S_M1_vrr = PAZ*I_ESP_Pz_S_M1_vrr-PRZ*I_ESP_Pz_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         ************************************************************/
        Double I_ESP_F3x_S_M1_vrr = PAX*I_ESP_D2x_S_M1_vrr-PRX*I_ESP_D2x_S_M2_vrr+2*oned2z*I_ESP_Px_S_M1_vrr-2*oned2z*I_ESP_Px_S_M2_vrr;
        Double I_ESP_F2xy_S_M1_vrr = PAY*I_ESP_D2x_S_M1_vrr-PRY*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_F2xz_S_M1_vrr = PAZ*I_ESP_D2x_S_M1_vrr-PRZ*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_Fx2y_S_M1_vrr = PAX*I_ESP_D2y_S_M1_vrr-PRX*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_Fx2z_S_M1_vrr = PAX*I_ESP_D2z_S_M1_vrr-PRX*I_ESP_D2z_S_M2_vrr;
        Double I_ESP_F3y_S_M1_vrr = PAY*I_ESP_D2y_S_M1_vrr-PRY*I_ESP_D2y_S_M2_vrr+2*oned2z*I_ESP_Py_S_M1_vrr-2*oned2z*I_ESP_Py_S_M2_vrr;
        Double I_ESP_F2yz_S_M1_vrr = PAZ*I_ESP_D2y_S_M1_vrr-PRZ*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_F3z_S_M1_vrr = PAZ*I_ESP_D2z_S_M1_vrr-PRZ*I_ESP_D2z_S_M2_vrr+2*oned2z*I_ESP_Pz_S_M1_vrr-2*oned2z*I_ESP_Pz_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_F_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_D_S_M2
         ************************************************************/
        Double I_ESP_G4x_S_M1_vrr = PAX*I_ESP_F3x_S_M1_vrr-PRX*I_ESP_F3x_S_M2_vrr+3*oned2z*I_ESP_D2x_S_M1_vrr-3*oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_G3xy_S_M1_vrr = PAY*I_ESP_F3x_S_M1_vrr-PRY*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_G3xz_S_M1_vrr = PAZ*I_ESP_F3x_S_M1_vrr-PRZ*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_G2x2y_S_M1_vrr = PAY*I_ESP_F2xy_S_M1_vrr-PRY*I_ESP_F2xy_S_M2_vrr+oned2z*I_ESP_D2x_S_M1_vrr-oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_G2x2z_S_M1_vrr = PAZ*I_ESP_F2xz_S_M1_vrr-PRZ*I_ESP_F2xz_S_M2_vrr+oned2z*I_ESP_D2x_S_M1_vrr-oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_Gx3y_S_M1_vrr = PAX*I_ESP_F3y_S_M1_vrr-PRX*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_Gx3z_S_M1_vrr = PAX*I_ESP_F3z_S_M1_vrr-PRX*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_G4y_S_M1_vrr = PAY*I_ESP_F3y_S_M1_vrr-PRY*I_ESP_F3y_S_M2_vrr+3*oned2z*I_ESP_D2y_S_M1_vrr-3*oned2z*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_G3yz_S_M1_vrr = PAZ*I_ESP_F3y_S_M1_vrr-PRZ*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_G2y2z_S_M1_vrr = PAZ*I_ESP_F2yz_S_M1_vrr-PRZ*I_ESP_F2yz_S_M2_vrr+oned2z*I_ESP_D2y_S_M1_vrr-oned2z*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_Gy3z_S_M1_vrr = PAY*I_ESP_F3z_S_M1_vrr-PRY*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_G4z_S_M1_vrr = PAZ*I_ESP_F3z_S_M1_vrr-PRZ*I_ESP_F3z_S_M2_vrr+3*oned2z*I_ESP_D2z_S_M1_vrr-3*oned2z*I_ESP_D2z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 5 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M1
         * RHS shell quartet name: SQ_ESP_G_S_M2
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_F_S_M2
         ************************************************************/
        Double I_ESP_H5x_S_M1_vrr = PAX*I_ESP_G4x_S_M1_vrr-PRX*I_ESP_G4x_S_M2_vrr+4*oned2z*I_ESP_F3x_S_M1_vrr-4*oned2z*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_H4xy_S_M1_vrr = PAY*I_ESP_G4x_S_M1_vrr-PRY*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_H4xz_S_M1_vrr = PAZ*I_ESP_G4x_S_M1_vrr-PRZ*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_H3x2y_S_M1_vrr = PAY*I_ESP_G3xy_S_M1_vrr-PRY*I_ESP_G3xy_S_M2_vrr+oned2z*I_ESP_F3x_S_M1_vrr-oned2z*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_H3x2z_S_M1_vrr = PAZ*I_ESP_G3xz_S_M1_vrr-PRZ*I_ESP_G3xz_S_M2_vrr+oned2z*I_ESP_F3x_S_M1_vrr-oned2z*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_H2x3y_S_M1_vrr = PAX*I_ESP_Gx3y_S_M1_vrr-PRX*I_ESP_Gx3y_S_M2_vrr+oned2z*I_ESP_F3y_S_M1_vrr-oned2z*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_H2x2yz_S_M1_vrr = PAZ*I_ESP_G2x2y_S_M1_vrr-PRZ*I_ESP_G2x2y_S_M2_vrr;
        Double I_ESP_H2x3z_S_M1_vrr = PAX*I_ESP_Gx3z_S_M1_vrr-PRX*I_ESP_Gx3z_S_M2_vrr+oned2z*I_ESP_F3z_S_M1_vrr-oned2z*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_Hx4y_S_M1_vrr = PAX*I_ESP_G4y_S_M1_vrr-PRX*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_Hx4z_S_M1_vrr = PAX*I_ESP_G4z_S_M1_vrr-PRX*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_H5y_S_M1_vrr = PAY*I_ESP_G4y_S_M1_vrr-PRY*I_ESP_G4y_S_M2_vrr+4*oned2z*I_ESP_F3y_S_M1_vrr-4*oned2z*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_H4yz_S_M1_vrr = PAZ*I_ESP_G4y_S_M1_vrr-PRZ*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_H3y2z_S_M1_vrr = PAZ*I_ESP_G3yz_S_M1_vrr-PRZ*I_ESP_G3yz_S_M2_vrr+oned2z*I_ESP_F3y_S_M1_vrr-oned2z*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_H2y3z_S_M1_vrr = PAY*I_ESP_Gy3z_S_M1_vrr-PRY*I_ESP_Gy3z_S_M2_vrr+oned2z*I_ESP_F3z_S_M1_vrr-oned2z*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_Hy4z_S_M1_vrr = PAY*I_ESP_G4z_S_M1_vrr-PRY*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_H5z_S_M1_vrr = PAZ*I_ESP_G4z_S_M1_vrr-PRZ*I_ESP_G4z_S_M2_vrr+4*oned2z*I_ESP_F3z_S_M1_vrr-4*oned2z*I_ESP_F3z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_Px_S_vrr = PAX*I_ESP_S_S_vrr-PRX*I_ESP_S_S_M1_vrr;
        Double I_ESP_Py_S_vrr = PAY*I_ESP_S_S_vrr-PRY*I_ESP_S_S_M1_vrr;
        Double I_ESP_Pz_S_vrr = PAZ*I_ESP_S_S_vrr-PRZ*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_D2x_S_vrr = PAX*I_ESP_Px_S_vrr-PRX*I_ESP_Px_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_Dxy_S_vrr = PAY*I_ESP_Px_S_vrr-PRY*I_ESP_Px_S_M1_vrr;
        Double I_ESP_Dxz_S_vrr = PAZ*I_ESP_Px_S_vrr-PRZ*I_ESP_Px_S_M1_vrr;
        Double I_ESP_D2y_S_vrr = PAY*I_ESP_Py_S_vrr-PRY*I_ESP_Py_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_Dyz_S_vrr = PAZ*I_ESP_Py_S_vrr-PRZ*I_ESP_Py_S_M1_vrr;
        Double I_ESP_D2z_S_vrr = PAZ*I_ESP_Pz_S_vrr-PRZ*I_ESP_Pz_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         ************************************************************/
        Double I_ESP_F3x_S_vrr = PAX*I_ESP_D2x_S_vrr-PRX*I_ESP_D2x_S_M1_vrr+2*oned2z*I_ESP_Px_S_vrr-2*oned2z*I_ESP_Px_S_M1_vrr;
        Double I_ESP_F2xy_S_vrr = PAY*I_ESP_D2x_S_vrr-PRY*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_F2xz_S_vrr = PAZ*I_ESP_D2x_S_vrr-PRZ*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Fx2y_S_vrr = PAX*I_ESP_D2y_S_vrr-PRX*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Fxyz_S_vrr = PAZ*I_ESP_Dxy_S_vrr-PRZ*I_ESP_Dxy_S_M1_vrr;
        Double I_ESP_Fx2z_S_vrr = PAX*I_ESP_D2z_S_vrr-PRX*I_ESP_D2z_S_M1_vrr;
        Double I_ESP_F3y_S_vrr = PAY*I_ESP_D2y_S_vrr-PRY*I_ESP_D2y_S_M1_vrr+2*oned2z*I_ESP_Py_S_vrr-2*oned2z*I_ESP_Py_S_M1_vrr;
        Double I_ESP_F2yz_S_vrr = PAZ*I_ESP_D2y_S_vrr-PRZ*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Fy2z_S_vrr = PAY*I_ESP_D2z_S_vrr-PRY*I_ESP_D2z_S_M1_vrr;
        Double I_ESP_F3z_S_vrr = PAZ*I_ESP_D2z_S_vrr-PRZ*I_ESP_D2z_S_M1_vrr+2*oned2z*I_ESP_Pz_S_vrr-2*oned2z*I_ESP_Pz_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         ************************************************************/
        Double I_ESP_G4x_S_vrr = PAX*I_ESP_F3x_S_vrr-PRX*I_ESP_F3x_S_M1_vrr+3*oned2z*I_ESP_D2x_S_vrr-3*oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G3xy_S_vrr = PAY*I_ESP_F3x_S_vrr-PRY*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G3xz_S_vrr = PAZ*I_ESP_F3x_S_vrr-PRZ*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G2x2y_S_vrr = PAY*I_ESP_F2xy_S_vrr-PRY*I_ESP_F2xy_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G2xyz_S_vrr = PAZ*I_ESP_F2xy_S_vrr-PRZ*I_ESP_F2xy_S_M1_vrr;
        Double I_ESP_G2x2z_S_vrr = PAZ*I_ESP_F2xz_S_vrr-PRZ*I_ESP_F2xz_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Gx3y_S_vrr = PAX*I_ESP_F3y_S_vrr-PRX*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_Gx2yz_S_vrr = PAZ*I_ESP_Fx2y_S_vrr-PRZ*I_ESP_Fx2y_S_M1_vrr;
        Double I_ESP_Gxy2z_S_vrr = PAY*I_ESP_Fx2z_S_vrr-PRY*I_ESP_Fx2z_S_M1_vrr;
        Double I_ESP_Gx3z_S_vrr = PAX*I_ESP_F3z_S_vrr-PRX*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_G4y_S_vrr = PAY*I_ESP_F3y_S_vrr-PRY*I_ESP_F3y_S_M1_vrr+3*oned2z*I_ESP_D2y_S_vrr-3*oned2z*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_G3yz_S_vrr = PAZ*I_ESP_F3y_S_vrr-PRZ*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_G2y2z_S_vrr = PAZ*I_ESP_F2yz_S_vrr-PRZ*I_ESP_F2yz_S_M1_vrr+oned2z*I_ESP_D2y_S_vrr-oned2z*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Gy3z_S_vrr = PAY*I_ESP_F3z_S_vrr-PRY*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_G4z_S_vrr = PAZ*I_ESP_F3z_S_vrr-PRZ*I_ESP_F3z_S_M1_vrr+3*oned2z*I_ESP_D2z_S_vrr-3*oned2z*I_ESP_D2z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S
         * RHS shell quartet name: SQ_ESP_G_S_M1
         * RHS shell quartet name: SQ_ESP_F_S
         * RHS shell quartet name: SQ_ESP_F_S_M1
         ************************************************************/
        Double I_ESP_H5x_S_vrr = PAX*I_ESP_G4x_S_vrr-PRX*I_ESP_G4x_S_M1_vrr+4*oned2z*I_ESP_F3x_S_vrr-4*oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H4xy_S_vrr = PAY*I_ESP_G4x_S_vrr-PRY*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_H4xz_S_vrr = PAZ*I_ESP_G4x_S_vrr-PRZ*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_H3x2y_S_vrr = PAY*I_ESP_G3xy_S_vrr-PRY*I_ESP_G3xy_S_M1_vrr+oned2z*I_ESP_F3x_S_vrr-oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H3xyz_S_vrr = PAZ*I_ESP_G3xy_S_vrr-PRZ*I_ESP_G3xy_S_M1_vrr;
        Double I_ESP_H3x2z_S_vrr = PAZ*I_ESP_G3xz_S_vrr-PRZ*I_ESP_G3xz_S_M1_vrr+oned2z*I_ESP_F3x_S_vrr-oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H2x3y_S_vrr = PAX*I_ESP_Gx3y_S_vrr-PRX*I_ESP_Gx3y_S_M1_vrr+oned2z*I_ESP_F3y_S_vrr-oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H2x2yz_S_vrr = PAZ*I_ESP_G2x2y_S_vrr-PRZ*I_ESP_G2x2y_S_M1_vrr;
        Double I_ESP_H2xy2z_S_vrr = PAY*I_ESP_G2x2z_S_vrr-PRY*I_ESP_G2x2z_S_M1_vrr;
        Double I_ESP_H2x3z_S_vrr = PAX*I_ESP_Gx3z_S_vrr-PRX*I_ESP_Gx3z_S_M1_vrr+oned2z*I_ESP_F3z_S_vrr-oned2z*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_Hx4y_S_vrr = PAX*I_ESP_G4y_S_vrr-PRX*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_Hx3yz_S_vrr = PAZ*I_ESP_Gx3y_S_vrr-PRZ*I_ESP_Gx3y_S_M1_vrr;
        Double I_ESP_Hx2y2z_S_vrr = PAX*I_ESP_G2y2z_S_vrr-PRX*I_ESP_G2y2z_S_M1_vrr;
        Double I_ESP_Hxy3z_S_vrr = PAY*I_ESP_Gx3z_S_vrr-PRY*I_ESP_Gx3z_S_M1_vrr;
        Double I_ESP_Hx4z_S_vrr = PAX*I_ESP_G4z_S_vrr-PRX*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_H5y_S_vrr = PAY*I_ESP_G4y_S_vrr-PRY*I_ESP_G4y_S_M1_vrr+4*oned2z*I_ESP_F3y_S_vrr-4*oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H4yz_S_vrr = PAZ*I_ESP_G4y_S_vrr-PRZ*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_H3y2z_S_vrr = PAZ*I_ESP_G3yz_S_vrr-PRZ*I_ESP_G3yz_S_M1_vrr+oned2z*I_ESP_F3y_S_vrr-oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H2y3z_S_vrr = PAY*I_ESP_Gy3z_S_vrr-PRY*I_ESP_Gy3z_S_M1_vrr+oned2z*I_ESP_F3z_S_vrr-oned2z*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_Hy4z_S_vrr = PAY*I_ESP_G4z_S_vrr-PRY*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_H5z_S_vrr = PAZ*I_ESP_G4z_S_vrr-PRZ*I_ESP_G4z_S_M1_vrr+4*oned2z*I_ESP_F3z_S_vrr-4*oned2z*I_ESP_F3z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S
         * RHS shell quartet name: SQ_ESP_H_S_M1
         * RHS shell quartet name: SQ_ESP_G_S
         * RHS shell quartet name: SQ_ESP_G_S_M1
         ************************************************************/
        Double I_ESP_I6x_S_vrr = PAX*I_ESP_H5x_S_vrr-PRX*I_ESP_H5x_S_M1_vrr+5*oned2z*I_ESP_G4x_S_vrr-5*oned2z*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_I5xy_S_vrr = PAY*I_ESP_H5x_S_vrr-PRY*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_I5xz_S_vrr = PAZ*I_ESP_H5x_S_vrr-PRZ*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_I4x2y_S_vrr = PAY*I_ESP_H4xy_S_vrr-PRY*I_ESP_H4xy_S_M1_vrr+oned2z*I_ESP_G4x_S_vrr-oned2z*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_I4xyz_S_vrr = PAZ*I_ESP_H4xy_S_vrr-PRZ*I_ESP_H4xy_S_M1_vrr;
        Double I_ESP_I4x2z_S_vrr = PAZ*I_ESP_H4xz_S_vrr-PRZ*I_ESP_H4xz_S_M1_vrr+oned2z*I_ESP_G4x_S_vrr-oned2z*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_I3x3y_S_vrr = PAY*I_ESP_H3x2y_S_vrr-PRY*I_ESP_H3x2y_S_M1_vrr+2*oned2z*I_ESP_G3xy_S_vrr-2*oned2z*I_ESP_G3xy_S_M1_vrr;
        Double I_ESP_I3x2yz_S_vrr = PAZ*I_ESP_H3x2y_S_vrr-PRZ*I_ESP_H3x2y_S_M1_vrr;
        Double I_ESP_I3xy2z_S_vrr = PAY*I_ESP_H3x2z_S_vrr-PRY*I_ESP_H3x2z_S_M1_vrr;
        Double I_ESP_I3x3z_S_vrr = PAZ*I_ESP_H3x2z_S_vrr-PRZ*I_ESP_H3x2z_S_M1_vrr+2*oned2z*I_ESP_G3xz_S_vrr-2*oned2z*I_ESP_G3xz_S_M1_vrr;
        Double I_ESP_I2x4y_S_vrr = PAX*I_ESP_Hx4y_S_vrr-PRX*I_ESP_Hx4y_S_M1_vrr+oned2z*I_ESP_G4y_S_vrr-oned2z*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_I2x3yz_S_vrr = PAZ*I_ESP_H2x3y_S_vrr-PRZ*I_ESP_H2x3y_S_M1_vrr;
        Double I_ESP_I2x2y2z_S_vrr = PAZ*I_ESP_H2x2yz_S_vrr-PRZ*I_ESP_H2x2yz_S_M1_vrr+oned2z*I_ESP_G2x2y_S_vrr-oned2z*I_ESP_G2x2y_S_M1_vrr;
        Double I_ESP_I2xy3z_S_vrr = PAY*I_ESP_H2x3z_S_vrr-PRY*I_ESP_H2x3z_S_M1_vrr;
        Double I_ESP_I2x4z_S_vrr = PAX*I_ESP_Hx4z_S_vrr-PRX*I_ESP_Hx4z_S_M1_vrr+oned2z*I_ESP_G4z_S_vrr-oned2z*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_Ix5y_S_vrr = PAX*I_ESP_H5y_S_vrr-PRX*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_Ix4yz_S_vrr = PAZ*I_ESP_Hx4y_S_vrr-PRZ*I_ESP_Hx4y_S_M1_vrr;
        Double I_ESP_Ix3y2z_S_vrr = PAX*I_ESP_H3y2z_S_vrr-PRX*I_ESP_H3y2z_S_M1_vrr;
        Double I_ESP_Ix2y3z_S_vrr = PAX*I_ESP_H2y3z_S_vrr-PRX*I_ESP_H2y3z_S_M1_vrr;
        Double I_ESP_Ixy4z_S_vrr = PAY*I_ESP_Hx4z_S_vrr-PRY*I_ESP_Hx4z_S_M1_vrr;
        Double I_ESP_Ix5z_S_vrr = PAX*I_ESP_H5z_S_vrr-PRX*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_I6y_S_vrr = PAY*I_ESP_H5y_S_vrr-PRY*I_ESP_H5y_S_M1_vrr+5*oned2z*I_ESP_G4y_S_vrr-5*oned2z*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_I5yz_S_vrr = PAZ*I_ESP_H5y_S_vrr-PRZ*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_I4y2z_S_vrr = PAZ*I_ESP_H4yz_S_vrr-PRZ*I_ESP_H4yz_S_M1_vrr+oned2z*I_ESP_G4y_S_vrr-oned2z*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_I3y3z_S_vrr = PAZ*I_ESP_H3y2z_S_vrr-PRZ*I_ESP_H3y2z_S_M1_vrr+2*oned2z*I_ESP_G3yz_S_vrr-2*oned2z*I_ESP_G3yz_S_M1_vrr;
        Double I_ESP_I2y4z_S_vrr = PAY*I_ESP_Hy4z_S_vrr-PRY*I_ESP_Hy4z_S_M1_vrr+oned2z*I_ESP_G4z_S_vrr-oned2z*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_Iy5z_S_vrr = PAY*I_ESP_H5z_S_vrr-PRY*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_I6z_S_vrr = PAZ*I_ESP_H5z_S_vrr-PRZ*I_ESP_H5z_S_M1_vrr+5*oned2z*I_ESP_G4z_S_vrr-5*oned2z*I_ESP_G4z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_I_S_aa_coefs = alpha*alpha;
        I_ESP_I6x_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I6x_S_vrr;
        I_ESP_I5xy_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I5xy_S_vrr;
        I_ESP_I5xz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I5xz_S_vrr;
        I_ESP_I4x2y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4x2y_S_vrr;
        I_ESP_I4xyz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4xyz_S_vrr;
        I_ESP_I4x2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4x2z_S_vrr;
        I_ESP_I3x3y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3x3y_S_vrr;
        I_ESP_I3x2yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3x2yz_S_vrr;
        I_ESP_I3xy2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3xy2z_S_vrr;
        I_ESP_I3x3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3x3z_S_vrr;
        I_ESP_I2x4y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x4y_S_vrr;
        I_ESP_I2x3yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x3yz_S_vrr;
        I_ESP_I2x2y2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x2y2z_S_vrr;
        I_ESP_I2xy3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2xy3z_S_vrr;
        I_ESP_I2x4z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x4z_S_vrr;
        I_ESP_Ix5y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix5y_S_vrr;
        I_ESP_Ix4yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix4yz_S_vrr;
        I_ESP_Ix3y2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix3y2z_S_vrr;
        I_ESP_Ix2y3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix2y3z_S_vrr;
        I_ESP_Ixy4z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ixy4z_S_vrr;
        I_ESP_Ix5z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix5z_S_vrr;
        I_ESP_I6y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I6y_S_vrr;
        I_ESP_I5yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I5yz_S_vrr;
        I_ESP_I4y2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4y2z_S_vrr;
        I_ESP_I3y3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3y3z_S_vrr;
        I_ESP_I2y4z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2y4z_S_vrr;
        I_ESP_Iy5z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Iy5z_S_vrr;
        I_ESP_I6z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I6z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_H_S_aa_coefs = alpha*alpha;
        I_ESP_H5x_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H5x_S_vrr;
        I_ESP_H4xy_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H4xy_S_vrr;
        I_ESP_H4xz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H4xz_S_vrr;
        I_ESP_H3x2y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3x2y_S_vrr;
        I_ESP_H3xyz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3xyz_S_vrr;
        I_ESP_H3x2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3x2z_S_vrr;
        I_ESP_H2x3y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2x3y_S_vrr;
        I_ESP_H2x2yz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2x2yz_S_vrr;
        I_ESP_H2xy2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2xy2z_S_vrr;
        I_ESP_H2x3z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2x3z_S_vrr;
        I_ESP_Hx4y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx4y_S_vrr;
        I_ESP_Hx3yz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx3yz_S_vrr;
        I_ESP_Hx2y2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx2y2z_S_vrr;
        I_ESP_Hxy3z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hxy3z_S_vrr;
        I_ESP_Hx4z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx4z_S_vrr;
        I_ESP_H5y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H5y_S_vrr;
        I_ESP_H4yz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H4yz_S_vrr;
        I_ESP_H3y2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3y2z_S_vrr;
        I_ESP_H2y3z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2y3z_S_vrr;
        I_ESP_Hy4z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hy4z_S_vrr;
        I_ESP_H5z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H5z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_G_S_a_coefs = alpha;
        I_ESP_G4x_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G4x_S_vrr;
        I_ESP_G3xy_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G3xy_S_vrr;
        I_ESP_G3xz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G3xz_S_vrr;
        I_ESP_G2x2y_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2x2y_S_vrr;
        I_ESP_G2xyz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2xyz_S_vrr;
        I_ESP_G2x2z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2x2z_S_vrr;
        I_ESP_Gx3y_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gx3y_S_vrr;
        I_ESP_Gx2yz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gx2yz_S_vrr;
        I_ESP_Gxy2z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gxy2z_S_vrr;
        I_ESP_Gx3z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gx3z_S_vrr;
        I_ESP_G4y_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G4y_S_vrr;
        I_ESP_G3yz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G3yz_S_vrr;
        I_ESP_G2y2z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2y2z_S_vrr;
        I_ESP_Gy3z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gy3z_S_vrr;
        I_ESP_G4z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G4z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_F_S_a_coefs = alpha;
        I_ESP_F3x_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F3x_S_vrr;
        I_ESP_F2xy_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F2xy_S_vrr;
        I_ESP_F2xz_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F2xz_S_vrr;
        I_ESP_Fx2y_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fx2y_S_vrr;
        I_ESP_Fxyz_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fxyz_S_vrr;
        I_ESP_Fx2z_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fx2z_S_vrr;
        I_ESP_F3y_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F3y_S_vrr;
        I_ESP_F2yz_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F2yz_S_vrr;
        I_ESP_Fy2z_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fy2z_S_vrr;
        I_ESP_F3z_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F3z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_D2x_S += I_ESP_D2x_S_vrr;
        I_ESP_Dxy_S += I_ESP_Dxy_S_vrr;
        I_ESP_Dxz_S += I_ESP_Dxz_S_vrr;
        I_ESP_D2y_S += I_ESP_D2y_S_vrr;
        I_ESP_Dyz_S += I_ESP_Dyz_S_vrr;
        I_ESP_D2z_S += I_ESP_D2z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_Px_S += I_ESP_Px_S_vrr;
        I_ESP_Py_S += I_ESP_Py_S_vrr;
        I_ESP_Pz_S += I_ESP_Pz_S_vrr;
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
     * shell quartet name: SQ_ESP_P_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_D_S
     * RHS shell quartet name: SQ_ESP_P_S
     ************************************************************/
    Double I_ESP_Px_Px = I_ESP_D2x_S+ABX*I_ESP_Px_S;
    Double I_ESP_Py_Px = I_ESP_Dxy_S+ABX*I_ESP_Py_S;
    Double I_ESP_Pz_Px = I_ESP_Dxz_S+ABX*I_ESP_Pz_S;
    Double I_ESP_Px_Py = I_ESP_Dxy_S+ABY*I_ESP_Px_S;
    Double I_ESP_Py_Py = I_ESP_D2y_S+ABY*I_ESP_Py_S;
    Double I_ESP_Pz_Py = I_ESP_Dyz_S+ABY*I_ESP_Pz_S;
    Double I_ESP_Px_Pz = I_ESP_Dxz_S+ABZ*I_ESP_Px_S;
    Double I_ESP_Py_Pz = I_ESP_Dyz_S+ABZ*I_ESP_Py_S;
    Double I_ESP_Pz_Pz = I_ESP_D2z_S+ABZ*I_ESP_Pz_S;

    /************************************************************
     * shell quartet name: SQ_ESP_F_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_S_a
     * RHS shell quartet name: SQ_ESP_F_S_a
     ************************************************************/
    Double I_ESP_F3x_Px_a = I_ESP_G4x_S_a+ABX*I_ESP_F3x_S_a;
    Double I_ESP_F2xy_Px_a = I_ESP_G3xy_S_a+ABX*I_ESP_F2xy_S_a;
    Double I_ESP_F2xz_Px_a = I_ESP_G3xz_S_a+ABX*I_ESP_F2xz_S_a;
    Double I_ESP_Fx2y_Px_a = I_ESP_G2x2y_S_a+ABX*I_ESP_Fx2y_S_a;
    Double I_ESP_Fxyz_Px_a = I_ESP_G2xyz_S_a+ABX*I_ESP_Fxyz_S_a;
    Double I_ESP_Fx2z_Px_a = I_ESP_G2x2z_S_a+ABX*I_ESP_Fx2z_S_a;
    Double I_ESP_F3y_Px_a = I_ESP_Gx3y_S_a+ABX*I_ESP_F3y_S_a;
    Double I_ESP_F2yz_Px_a = I_ESP_Gx2yz_S_a+ABX*I_ESP_F2yz_S_a;
    Double I_ESP_Fy2z_Px_a = I_ESP_Gxy2z_S_a+ABX*I_ESP_Fy2z_S_a;
    Double I_ESP_F3z_Px_a = I_ESP_Gx3z_S_a+ABX*I_ESP_F3z_S_a;
    Double I_ESP_F3x_Py_a = I_ESP_G3xy_S_a+ABY*I_ESP_F3x_S_a;
    Double I_ESP_F2xy_Py_a = I_ESP_G2x2y_S_a+ABY*I_ESP_F2xy_S_a;
    Double I_ESP_F2xz_Py_a = I_ESP_G2xyz_S_a+ABY*I_ESP_F2xz_S_a;
    Double I_ESP_Fx2y_Py_a = I_ESP_Gx3y_S_a+ABY*I_ESP_Fx2y_S_a;
    Double I_ESP_Fxyz_Py_a = I_ESP_Gx2yz_S_a+ABY*I_ESP_Fxyz_S_a;
    Double I_ESP_Fx2z_Py_a = I_ESP_Gxy2z_S_a+ABY*I_ESP_Fx2z_S_a;
    Double I_ESP_F3y_Py_a = I_ESP_G4y_S_a+ABY*I_ESP_F3y_S_a;
    Double I_ESP_F2yz_Py_a = I_ESP_G3yz_S_a+ABY*I_ESP_F2yz_S_a;
    Double I_ESP_Fy2z_Py_a = I_ESP_G2y2z_S_a+ABY*I_ESP_Fy2z_S_a;
    Double I_ESP_F3z_Py_a = I_ESP_Gy3z_S_a+ABY*I_ESP_F3z_S_a;
    Double I_ESP_F3x_Pz_a = I_ESP_G3xz_S_a+ABZ*I_ESP_F3x_S_a;
    Double I_ESP_F2xy_Pz_a = I_ESP_G2xyz_S_a+ABZ*I_ESP_F2xy_S_a;
    Double I_ESP_F2xz_Pz_a = I_ESP_G2x2z_S_a+ABZ*I_ESP_F2xz_S_a;
    Double I_ESP_Fx2y_Pz_a = I_ESP_Gx2yz_S_a+ABZ*I_ESP_Fx2y_S_a;
    Double I_ESP_Fxyz_Pz_a = I_ESP_Gxy2z_S_a+ABZ*I_ESP_Fxyz_S_a;
    Double I_ESP_Fx2z_Pz_a = I_ESP_Gx3z_S_a+ABZ*I_ESP_Fx2z_S_a;
    Double I_ESP_F3y_Pz_a = I_ESP_G3yz_S_a+ABZ*I_ESP_F3y_S_a;
    Double I_ESP_F2yz_Pz_a = I_ESP_G2y2z_S_a+ABZ*I_ESP_F2yz_S_a;
    Double I_ESP_Fy2z_Pz_a = I_ESP_Gy3z_S_a+ABZ*I_ESP_Fy2z_S_a;
    Double I_ESP_F3z_Pz_a = I_ESP_G4z_S_a+ABZ*I_ESP_F3z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S_aa
     * RHS shell quartet name: SQ_ESP_H_S_aa
     ************************************************************/
    Double I_ESP_H5x_Px_aa = I_ESP_I6x_S_aa+ABX*I_ESP_H5x_S_aa;
    Double I_ESP_H4xy_Px_aa = I_ESP_I5xy_S_aa+ABX*I_ESP_H4xy_S_aa;
    Double I_ESP_H4xz_Px_aa = I_ESP_I5xz_S_aa+ABX*I_ESP_H4xz_S_aa;
    Double I_ESP_H3x2y_Px_aa = I_ESP_I4x2y_S_aa+ABX*I_ESP_H3x2y_S_aa;
    Double I_ESP_H3xyz_Px_aa = I_ESP_I4xyz_S_aa+ABX*I_ESP_H3xyz_S_aa;
    Double I_ESP_H3x2z_Px_aa = I_ESP_I4x2z_S_aa+ABX*I_ESP_H3x2z_S_aa;
    Double I_ESP_H2x3y_Px_aa = I_ESP_I3x3y_S_aa+ABX*I_ESP_H2x3y_S_aa;
    Double I_ESP_H2x2yz_Px_aa = I_ESP_I3x2yz_S_aa+ABX*I_ESP_H2x2yz_S_aa;
    Double I_ESP_H2xy2z_Px_aa = I_ESP_I3xy2z_S_aa+ABX*I_ESP_H2xy2z_S_aa;
    Double I_ESP_H2x3z_Px_aa = I_ESP_I3x3z_S_aa+ABX*I_ESP_H2x3z_S_aa;
    Double I_ESP_Hx4y_Px_aa = I_ESP_I2x4y_S_aa+ABX*I_ESP_Hx4y_S_aa;
    Double I_ESP_Hx3yz_Px_aa = I_ESP_I2x3yz_S_aa+ABX*I_ESP_Hx3yz_S_aa;
    Double I_ESP_Hx2y2z_Px_aa = I_ESP_I2x2y2z_S_aa+ABX*I_ESP_Hx2y2z_S_aa;
    Double I_ESP_Hxy3z_Px_aa = I_ESP_I2xy3z_S_aa+ABX*I_ESP_Hxy3z_S_aa;
    Double I_ESP_Hx4z_Px_aa = I_ESP_I2x4z_S_aa+ABX*I_ESP_Hx4z_S_aa;
    Double I_ESP_H5y_Px_aa = I_ESP_Ix5y_S_aa+ABX*I_ESP_H5y_S_aa;
    Double I_ESP_H4yz_Px_aa = I_ESP_Ix4yz_S_aa+ABX*I_ESP_H4yz_S_aa;
    Double I_ESP_H3y2z_Px_aa = I_ESP_Ix3y2z_S_aa+ABX*I_ESP_H3y2z_S_aa;
    Double I_ESP_H2y3z_Px_aa = I_ESP_Ix2y3z_S_aa+ABX*I_ESP_H2y3z_S_aa;
    Double I_ESP_Hy4z_Px_aa = I_ESP_Ixy4z_S_aa+ABX*I_ESP_Hy4z_S_aa;
    Double I_ESP_H5z_Px_aa = I_ESP_Ix5z_S_aa+ABX*I_ESP_H5z_S_aa;
    Double I_ESP_H5x_Py_aa = I_ESP_I5xy_S_aa+ABY*I_ESP_H5x_S_aa;
    Double I_ESP_H4xy_Py_aa = I_ESP_I4x2y_S_aa+ABY*I_ESP_H4xy_S_aa;
    Double I_ESP_H4xz_Py_aa = I_ESP_I4xyz_S_aa+ABY*I_ESP_H4xz_S_aa;
    Double I_ESP_H3x2y_Py_aa = I_ESP_I3x3y_S_aa+ABY*I_ESP_H3x2y_S_aa;
    Double I_ESP_H3xyz_Py_aa = I_ESP_I3x2yz_S_aa+ABY*I_ESP_H3xyz_S_aa;
    Double I_ESP_H3x2z_Py_aa = I_ESP_I3xy2z_S_aa+ABY*I_ESP_H3x2z_S_aa;
    Double I_ESP_H2x3y_Py_aa = I_ESP_I2x4y_S_aa+ABY*I_ESP_H2x3y_S_aa;
    Double I_ESP_H2x2yz_Py_aa = I_ESP_I2x3yz_S_aa+ABY*I_ESP_H2x2yz_S_aa;
    Double I_ESP_H2xy2z_Py_aa = I_ESP_I2x2y2z_S_aa+ABY*I_ESP_H2xy2z_S_aa;
    Double I_ESP_H2x3z_Py_aa = I_ESP_I2xy3z_S_aa+ABY*I_ESP_H2x3z_S_aa;
    Double I_ESP_Hx4y_Py_aa = I_ESP_Ix5y_S_aa+ABY*I_ESP_Hx4y_S_aa;
    Double I_ESP_Hx3yz_Py_aa = I_ESP_Ix4yz_S_aa+ABY*I_ESP_Hx3yz_S_aa;
    Double I_ESP_Hx2y2z_Py_aa = I_ESP_Ix3y2z_S_aa+ABY*I_ESP_Hx2y2z_S_aa;
    Double I_ESP_Hxy3z_Py_aa = I_ESP_Ix2y3z_S_aa+ABY*I_ESP_Hxy3z_S_aa;
    Double I_ESP_Hx4z_Py_aa = I_ESP_Ixy4z_S_aa+ABY*I_ESP_Hx4z_S_aa;
    Double I_ESP_H5y_Py_aa = I_ESP_I6y_S_aa+ABY*I_ESP_H5y_S_aa;
    Double I_ESP_H4yz_Py_aa = I_ESP_I5yz_S_aa+ABY*I_ESP_H4yz_S_aa;
    Double I_ESP_H3y2z_Py_aa = I_ESP_I4y2z_S_aa+ABY*I_ESP_H3y2z_S_aa;
    Double I_ESP_H2y3z_Py_aa = I_ESP_I3y3z_S_aa+ABY*I_ESP_H2y3z_S_aa;
    Double I_ESP_Hy4z_Py_aa = I_ESP_I2y4z_S_aa+ABY*I_ESP_Hy4z_S_aa;
    Double I_ESP_H5z_Py_aa = I_ESP_Iy5z_S_aa+ABY*I_ESP_H5z_S_aa;
    Double I_ESP_H5x_Pz_aa = I_ESP_I5xz_S_aa+ABZ*I_ESP_H5x_S_aa;
    Double I_ESP_H4xy_Pz_aa = I_ESP_I4xyz_S_aa+ABZ*I_ESP_H4xy_S_aa;
    Double I_ESP_H4xz_Pz_aa = I_ESP_I4x2z_S_aa+ABZ*I_ESP_H4xz_S_aa;
    Double I_ESP_H3x2y_Pz_aa = I_ESP_I3x2yz_S_aa+ABZ*I_ESP_H3x2y_S_aa;
    Double I_ESP_H3xyz_Pz_aa = I_ESP_I3xy2z_S_aa+ABZ*I_ESP_H3xyz_S_aa;
    Double I_ESP_H3x2z_Pz_aa = I_ESP_I3x3z_S_aa+ABZ*I_ESP_H3x2z_S_aa;
    Double I_ESP_H2x3y_Pz_aa = I_ESP_I2x3yz_S_aa+ABZ*I_ESP_H2x3y_S_aa;
    Double I_ESP_H2x2yz_Pz_aa = I_ESP_I2x2y2z_S_aa+ABZ*I_ESP_H2x2yz_S_aa;
    Double I_ESP_H2xy2z_Pz_aa = I_ESP_I2xy3z_S_aa+ABZ*I_ESP_H2xy2z_S_aa;
    Double I_ESP_H2x3z_Pz_aa = I_ESP_I2x4z_S_aa+ABZ*I_ESP_H2x3z_S_aa;
    Double I_ESP_Hx4y_Pz_aa = I_ESP_Ix4yz_S_aa+ABZ*I_ESP_Hx4y_S_aa;
    Double I_ESP_Hx3yz_Pz_aa = I_ESP_Ix3y2z_S_aa+ABZ*I_ESP_Hx3yz_S_aa;
    Double I_ESP_Hx2y2z_Pz_aa = I_ESP_Ix2y3z_S_aa+ABZ*I_ESP_Hx2y2z_S_aa;
    Double I_ESP_Hxy3z_Pz_aa = I_ESP_Ixy4z_S_aa+ABZ*I_ESP_Hxy3z_S_aa;
    Double I_ESP_Hx4z_Pz_aa = I_ESP_Ix5z_S_aa+ABZ*I_ESP_Hx4z_S_aa;
    Double I_ESP_H5y_Pz_aa = I_ESP_I5yz_S_aa+ABZ*I_ESP_H5y_S_aa;
    Double I_ESP_H4yz_Pz_aa = I_ESP_I4y2z_S_aa+ABZ*I_ESP_H4yz_S_aa;
    Double I_ESP_H3y2z_Pz_aa = I_ESP_I3y3z_S_aa+ABZ*I_ESP_H3y2z_S_aa;
    Double I_ESP_H2y3z_Pz_aa = I_ESP_I2y4z_S_aa+ABZ*I_ESP_H2y3z_S_aa;
    Double I_ESP_Hy4z_Pz_aa = I_ESP_Iy5z_S_aa+ABZ*I_ESP_Hy4z_S_aa;
    Double I_ESP_H5z_Pz_aa = I_ESP_I6z_S_aa+ABZ*I_ESP_H5z_S_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_F_P_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_P_aa
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_P_P
     ************************************************************/
    abcd[iGrid*180+0] = 4.0E0*I_ESP_H5x_Px_aa-2.0E0*3*I_ESP_F3x_Px_a-2.0E0*4*I_ESP_F3x_Px_a+3*2*I_ESP_Px_Px;
    abcd[iGrid*180+1] = 4.0E0*I_ESP_H4xy_Px_aa-2.0E0*2*I_ESP_F2xy_Px_a-2.0E0*3*I_ESP_F2xy_Px_a+2*1*I_ESP_Py_Px;
    abcd[iGrid*180+2] = 4.0E0*I_ESP_H4xz_Px_aa-2.0E0*2*I_ESP_F2xz_Px_a-2.0E0*3*I_ESP_F2xz_Px_a+2*1*I_ESP_Pz_Px;
    abcd[iGrid*180+3] = 4.0E0*I_ESP_H3x2y_Px_aa-2.0E0*1*I_ESP_Fx2y_Px_a-2.0E0*2*I_ESP_Fx2y_Px_a;
    abcd[iGrid*180+4] = 4.0E0*I_ESP_H3xyz_Px_aa-2.0E0*1*I_ESP_Fxyz_Px_a-2.0E0*2*I_ESP_Fxyz_Px_a;
    abcd[iGrid*180+5] = 4.0E0*I_ESP_H3x2z_Px_aa-2.0E0*1*I_ESP_Fx2z_Px_a-2.0E0*2*I_ESP_Fx2z_Px_a;
    abcd[iGrid*180+6] = 4.0E0*I_ESP_H2x3y_Px_aa-2.0E0*1*I_ESP_F3y_Px_a;
    abcd[iGrid*180+7] = 4.0E0*I_ESP_H2x2yz_Px_aa-2.0E0*1*I_ESP_F2yz_Px_a;
    abcd[iGrid*180+8] = 4.0E0*I_ESP_H2xy2z_Px_aa-2.0E0*1*I_ESP_Fy2z_Px_a;
    abcd[iGrid*180+9] = 4.0E0*I_ESP_H2x3z_Px_aa-2.0E0*1*I_ESP_F3z_Px_a;
    abcd[iGrid*180+10] = 4.0E0*I_ESP_H5x_Py_aa-2.0E0*3*I_ESP_F3x_Py_a-2.0E0*4*I_ESP_F3x_Py_a+3*2*I_ESP_Px_Py;
    abcd[iGrid*180+11] = 4.0E0*I_ESP_H4xy_Py_aa-2.0E0*2*I_ESP_F2xy_Py_a-2.0E0*3*I_ESP_F2xy_Py_a+2*1*I_ESP_Py_Py;
    abcd[iGrid*180+12] = 4.0E0*I_ESP_H4xz_Py_aa-2.0E0*2*I_ESP_F2xz_Py_a-2.0E0*3*I_ESP_F2xz_Py_a+2*1*I_ESP_Pz_Py;
    abcd[iGrid*180+13] = 4.0E0*I_ESP_H3x2y_Py_aa-2.0E0*1*I_ESP_Fx2y_Py_a-2.0E0*2*I_ESP_Fx2y_Py_a;
    abcd[iGrid*180+14] = 4.0E0*I_ESP_H3xyz_Py_aa-2.0E0*1*I_ESP_Fxyz_Py_a-2.0E0*2*I_ESP_Fxyz_Py_a;
    abcd[iGrid*180+15] = 4.0E0*I_ESP_H3x2z_Py_aa-2.0E0*1*I_ESP_Fx2z_Py_a-2.0E0*2*I_ESP_Fx2z_Py_a;
    abcd[iGrid*180+16] = 4.0E0*I_ESP_H2x3y_Py_aa-2.0E0*1*I_ESP_F3y_Py_a;
    abcd[iGrid*180+17] = 4.0E0*I_ESP_H2x2yz_Py_aa-2.0E0*1*I_ESP_F2yz_Py_a;
    abcd[iGrid*180+18] = 4.0E0*I_ESP_H2xy2z_Py_aa-2.0E0*1*I_ESP_Fy2z_Py_a;
    abcd[iGrid*180+19] = 4.0E0*I_ESP_H2x3z_Py_aa-2.0E0*1*I_ESP_F3z_Py_a;
    abcd[iGrid*180+20] = 4.0E0*I_ESP_H5x_Pz_aa-2.0E0*3*I_ESP_F3x_Pz_a-2.0E0*4*I_ESP_F3x_Pz_a+3*2*I_ESP_Px_Pz;
    abcd[iGrid*180+21] = 4.0E0*I_ESP_H4xy_Pz_aa-2.0E0*2*I_ESP_F2xy_Pz_a-2.0E0*3*I_ESP_F2xy_Pz_a+2*1*I_ESP_Py_Pz;
    abcd[iGrid*180+22] = 4.0E0*I_ESP_H4xz_Pz_aa-2.0E0*2*I_ESP_F2xz_Pz_a-2.0E0*3*I_ESP_F2xz_Pz_a+2*1*I_ESP_Pz_Pz;
    abcd[iGrid*180+23] = 4.0E0*I_ESP_H3x2y_Pz_aa-2.0E0*1*I_ESP_Fx2y_Pz_a-2.0E0*2*I_ESP_Fx2y_Pz_a;
    abcd[iGrid*180+24] = 4.0E0*I_ESP_H3xyz_Pz_aa-2.0E0*1*I_ESP_Fxyz_Pz_a-2.0E0*2*I_ESP_Fxyz_Pz_a;
    abcd[iGrid*180+25] = 4.0E0*I_ESP_H3x2z_Pz_aa-2.0E0*1*I_ESP_Fx2z_Pz_a-2.0E0*2*I_ESP_Fx2z_Pz_a;
    abcd[iGrid*180+26] = 4.0E0*I_ESP_H2x3y_Pz_aa-2.0E0*1*I_ESP_F3y_Pz_a;
    abcd[iGrid*180+27] = 4.0E0*I_ESP_H2x2yz_Pz_aa-2.0E0*1*I_ESP_F2yz_Pz_a;
    abcd[iGrid*180+28] = 4.0E0*I_ESP_H2xy2z_Pz_aa-2.0E0*1*I_ESP_Fy2z_Pz_a;
    abcd[iGrid*180+29] = 4.0E0*I_ESP_H2x3z_Pz_aa-2.0E0*1*I_ESP_F3z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_P_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_P_aa
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_P_P
     ************************************************************/
    abcd[iGrid*180+30] = 4.0E0*I_ESP_H4xy_Px_aa-2.0E0*3*I_ESP_F2xy_Px_a;
    abcd[iGrid*180+31] = 4.0E0*I_ESP_H3x2y_Px_aa-2.0E0*1*I_ESP_F3x_Px_a-2.0E0*2*I_ESP_Fx2y_Px_a+2*1*I_ESP_Px_Px;
    abcd[iGrid*180+32] = 4.0E0*I_ESP_H3xyz_Px_aa-2.0E0*2*I_ESP_Fxyz_Px_a;
    abcd[iGrid*180+33] = 4.0E0*I_ESP_H2x3y_Px_aa-2.0E0*2*I_ESP_F2xy_Px_a-2.0E0*1*I_ESP_F3y_Px_a+2*I_ESP_Py_Px;
    abcd[iGrid*180+34] = 4.0E0*I_ESP_H2x2yz_Px_aa-2.0E0*1*I_ESP_F2xz_Px_a-2.0E0*1*I_ESP_F2yz_Px_a+1*I_ESP_Pz_Px;
    abcd[iGrid*180+35] = 4.0E0*I_ESP_H2xy2z_Px_aa-2.0E0*1*I_ESP_Fy2z_Px_a;
    abcd[iGrid*180+36] = 4.0E0*I_ESP_Hx4y_Px_aa-2.0E0*3*I_ESP_Fx2y_Px_a;
    abcd[iGrid*180+37] = 4.0E0*I_ESP_Hx3yz_Px_aa-2.0E0*2*I_ESP_Fxyz_Px_a;
    abcd[iGrid*180+38] = 4.0E0*I_ESP_Hx2y2z_Px_aa-2.0E0*1*I_ESP_Fx2z_Px_a;
    abcd[iGrid*180+39] = 4.0E0*I_ESP_Hxy3z_Px_aa;
    abcd[iGrid*180+40] = 4.0E0*I_ESP_H4xy_Py_aa-2.0E0*3*I_ESP_F2xy_Py_a;
    abcd[iGrid*180+41] = 4.0E0*I_ESP_H3x2y_Py_aa-2.0E0*1*I_ESP_F3x_Py_a-2.0E0*2*I_ESP_Fx2y_Py_a+2*1*I_ESP_Px_Py;
    abcd[iGrid*180+42] = 4.0E0*I_ESP_H3xyz_Py_aa-2.0E0*2*I_ESP_Fxyz_Py_a;
    abcd[iGrid*180+43] = 4.0E0*I_ESP_H2x3y_Py_aa-2.0E0*2*I_ESP_F2xy_Py_a-2.0E0*1*I_ESP_F3y_Py_a+2*I_ESP_Py_Py;
    abcd[iGrid*180+44] = 4.0E0*I_ESP_H2x2yz_Py_aa-2.0E0*1*I_ESP_F2xz_Py_a-2.0E0*1*I_ESP_F2yz_Py_a+1*I_ESP_Pz_Py;
    abcd[iGrid*180+45] = 4.0E0*I_ESP_H2xy2z_Py_aa-2.0E0*1*I_ESP_Fy2z_Py_a;
    abcd[iGrid*180+46] = 4.0E0*I_ESP_Hx4y_Py_aa-2.0E0*3*I_ESP_Fx2y_Py_a;
    abcd[iGrid*180+47] = 4.0E0*I_ESP_Hx3yz_Py_aa-2.0E0*2*I_ESP_Fxyz_Py_a;
    abcd[iGrid*180+48] = 4.0E0*I_ESP_Hx2y2z_Py_aa-2.0E0*1*I_ESP_Fx2z_Py_a;
    abcd[iGrid*180+49] = 4.0E0*I_ESP_Hxy3z_Py_aa;
    abcd[iGrid*180+50] = 4.0E0*I_ESP_H4xy_Pz_aa-2.0E0*3*I_ESP_F2xy_Pz_a;
    abcd[iGrid*180+51] = 4.0E0*I_ESP_H3x2y_Pz_aa-2.0E0*1*I_ESP_F3x_Pz_a-2.0E0*2*I_ESP_Fx2y_Pz_a+2*1*I_ESP_Px_Pz;
    abcd[iGrid*180+52] = 4.0E0*I_ESP_H3xyz_Pz_aa-2.0E0*2*I_ESP_Fxyz_Pz_a;
    abcd[iGrid*180+53] = 4.0E0*I_ESP_H2x3y_Pz_aa-2.0E0*2*I_ESP_F2xy_Pz_a-2.0E0*1*I_ESP_F3y_Pz_a+2*I_ESP_Py_Pz;
    abcd[iGrid*180+54] = 4.0E0*I_ESP_H2x2yz_Pz_aa-2.0E0*1*I_ESP_F2xz_Pz_a-2.0E0*1*I_ESP_F2yz_Pz_a+1*I_ESP_Pz_Pz;
    abcd[iGrid*180+55] = 4.0E0*I_ESP_H2xy2z_Pz_aa-2.0E0*1*I_ESP_Fy2z_Pz_a;
    abcd[iGrid*180+56] = 4.0E0*I_ESP_Hx4y_Pz_aa-2.0E0*3*I_ESP_Fx2y_Pz_a;
    abcd[iGrid*180+57] = 4.0E0*I_ESP_Hx3yz_Pz_aa-2.0E0*2*I_ESP_Fxyz_Pz_a;
    abcd[iGrid*180+58] = 4.0E0*I_ESP_Hx2y2z_Pz_aa-2.0E0*1*I_ESP_Fx2z_Pz_a;
    abcd[iGrid*180+59] = 4.0E0*I_ESP_Hxy3z_Pz_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_F_P_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_P_aa
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_P_P
     ************************************************************/
    abcd[iGrid*180+60] = 4.0E0*I_ESP_H4xz_Px_aa-2.0E0*3*I_ESP_F2xz_Px_a;
    abcd[iGrid*180+61] = 4.0E0*I_ESP_H3xyz_Px_aa-2.0E0*2*I_ESP_Fxyz_Px_a;
    abcd[iGrid*180+62] = 4.0E0*I_ESP_H3x2z_Px_aa-2.0E0*1*I_ESP_F3x_Px_a-2.0E0*2*I_ESP_Fx2z_Px_a+2*1*I_ESP_Px_Px;
    abcd[iGrid*180+63] = 4.0E0*I_ESP_H2x2yz_Px_aa-2.0E0*1*I_ESP_F2yz_Px_a;
    abcd[iGrid*180+64] = 4.0E0*I_ESP_H2xy2z_Px_aa-2.0E0*1*I_ESP_F2xy_Px_a-2.0E0*1*I_ESP_Fy2z_Px_a+1*I_ESP_Py_Px;
    abcd[iGrid*180+65] = 4.0E0*I_ESP_H2x3z_Px_aa-2.0E0*2*I_ESP_F2xz_Px_a-2.0E0*1*I_ESP_F3z_Px_a+2*I_ESP_Pz_Px;
    abcd[iGrid*180+66] = 4.0E0*I_ESP_Hx3yz_Px_aa;
    abcd[iGrid*180+67] = 4.0E0*I_ESP_Hx2y2z_Px_aa-2.0E0*1*I_ESP_Fx2y_Px_a;
    abcd[iGrid*180+68] = 4.0E0*I_ESP_Hxy3z_Px_aa-2.0E0*2*I_ESP_Fxyz_Px_a;
    abcd[iGrid*180+69] = 4.0E0*I_ESP_Hx4z_Px_aa-2.0E0*3*I_ESP_Fx2z_Px_a;
    abcd[iGrid*180+70] = 4.0E0*I_ESP_H4xz_Py_aa-2.0E0*3*I_ESP_F2xz_Py_a;
    abcd[iGrid*180+71] = 4.0E0*I_ESP_H3xyz_Py_aa-2.0E0*2*I_ESP_Fxyz_Py_a;
    abcd[iGrid*180+72] = 4.0E0*I_ESP_H3x2z_Py_aa-2.0E0*1*I_ESP_F3x_Py_a-2.0E0*2*I_ESP_Fx2z_Py_a+2*1*I_ESP_Px_Py;
    abcd[iGrid*180+73] = 4.0E0*I_ESP_H2x2yz_Py_aa-2.0E0*1*I_ESP_F2yz_Py_a;
    abcd[iGrid*180+74] = 4.0E0*I_ESP_H2xy2z_Py_aa-2.0E0*1*I_ESP_F2xy_Py_a-2.0E0*1*I_ESP_Fy2z_Py_a+1*I_ESP_Py_Py;
    abcd[iGrid*180+75] = 4.0E0*I_ESP_H2x3z_Py_aa-2.0E0*2*I_ESP_F2xz_Py_a-2.0E0*1*I_ESP_F3z_Py_a+2*I_ESP_Pz_Py;
    abcd[iGrid*180+76] = 4.0E0*I_ESP_Hx3yz_Py_aa;
    abcd[iGrid*180+77] = 4.0E0*I_ESP_Hx2y2z_Py_aa-2.0E0*1*I_ESP_Fx2y_Py_a;
    abcd[iGrid*180+78] = 4.0E0*I_ESP_Hxy3z_Py_aa-2.0E0*2*I_ESP_Fxyz_Py_a;
    abcd[iGrid*180+79] = 4.0E0*I_ESP_Hx4z_Py_aa-2.0E0*3*I_ESP_Fx2z_Py_a;
    abcd[iGrid*180+80] = 4.0E0*I_ESP_H4xz_Pz_aa-2.0E0*3*I_ESP_F2xz_Pz_a;
    abcd[iGrid*180+81] = 4.0E0*I_ESP_H3xyz_Pz_aa-2.0E0*2*I_ESP_Fxyz_Pz_a;
    abcd[iGrid*180+82] = 4.0E0*I_ESP_H3x2z_Pz_aa-2.0E0*1*I_ESP_F3x_Pz_a-2.0E0*2*I_ESP_Fx2z_Pz_a+2*1*I_ESP_Px_Pz;
    abcd[iGrid*180+83] = 4.0E0*I_ESP_H2x2yz_Pz_aa-2.0E0*1*I_ESP_F2yz_Pz_a;
    abcd[iGrid*180+84] = 4.0E0*I_ESP_H2xy2z_Pz_aa-2.0E0*1*I_ESP_F2xy_Pz_a-2.0E0*1*I_ESP_Fy2z_Pz_a+1*I_ESP_Py_Pz;
    abcd[iGrid*180+85] = 4.0E0*I_ESP_H2x3z_Pz_aa-2.0E0*2*I_ESP_F2xz_Pz_a-2.0E0*1*I_ESP_F3z_Pz_a+2*I_ESP_Pz_Pz;
    abcd[iGrid*180+86] = 4.0E0*I_ESP_Hx3yz_Pz_aa;
    abcd[iGrid*180+87] = 4.0E0*I_ESP_Hx2y2z_Pz_aa-2.0E0*1*I_ESP_Fx2y_Pz_a;
    abcd[iGrid*180+88] = 4.0E0*I_ESP_Hxy3z_Pz_aa-2.0E0*2*I_ESP_Fxyz_Pz_a;
    abcd[iGrid*180+89] = 4.0E0*I_ESP_Hx4z_Pz_aa-2.0E0*3*I_ESP_Fx2z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_P_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_P_aa
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_P_P
     ************************************************************/
    abcd[iGrid*180+90] = 4.0E0*I_ESP_H3x2y_Px_aa-2.0E0*1*I_ESP_F3x_Px_a;
    abcd[iGrid*180+91] = 4.0E0*I_ESP_H2x3y_Px_aa-2.0E0*1*I_ESP_F2xy_Px_a-2.0E0*2*I_ESP_F2xy_Px_a;
    abcd[iGrid*180+92] = 4.0E0*I_ESP_H2x2yz_Px_aa-2.0E0*1*I_ESP_F2xz_Px_a;
    abcd[iGrid*180+93] = 4.0E0*I_ESP_Hx4y_Px_aa-2.0E0*2*I_ESP_Fx2y_Px_a-2.0E0*3*I_ESP_Fx2y_Px_a+2*1*I_ESP_Px_Px;
    abcd[iGrid*180+94] = 4.0E0*I_ESP_Hx3yz_Px_aa-2.0E0*1*I_ESP_Fxyz_Px_a-2.0E0*2*I_ESP_Fxyz_Px_a;
    abcd[iGrid*180+95] = 4.0E0*I_ESP_Hx2y2z_Px_aa-2.0E0*1*I_ESP_Fx2z_Px_a;
    abcd[iGrid*180+96] = 4.0E0*I_ESP_H5y_Px_aa-2.0E0*3*I_ESP_F3y_Px_a-2.0E0*4*I_ESP_F3y_Px_a+3*2*I_ESP_Py_Px;
    abcd[iGrid*180+97] = 4.0E0*I_ESP_H4yz_Px_aa-2.0E0*2*I_ESP_F2yz_Px_a-2.0E0*3*I_ESP_F2yz_Px_a+2*1*I_ESP_Pz_Px;
    abcd[iGrid*180+98] = 4.0E0*I_ESP_H3y2z_Px_aa-2.0E0*1*I_ESP_Fy2z_Px_a-2.0E0*2*I_ESP_Fy2z_Px_a;
    abcd[iGrid*180+99] = 4.0E0*I_ESP_H2y3z_Px_aa-2.0E0*1*I_ESP_F3z_Px_a;
    abcd[iGrid*180+100] = 4.0E0*I_ESP_H3x2y_Py_aa-2.0E0*1*I_ESP_F3x_Py_a;
    abcd[iGrid*180+101] = 4.0E0*I_ESP_H2x3y_Py_aa-2.0E0*1*I_ESP_F2xy_Py_a-2.0E0*2*I_ESP_F2xy_Py_a;
    abcd[iGrid*180+102] = 4.0E0*I_ESP_H2x2yz_Py_aa-2.0E0*1*I_ESP_F2xz_Py_a;
    abcd[iGrid*180+103] = 4.0E0*I_ESP_Hx4y_Py_aa-2.0E0*2*I_ESP_Fx2y_Py_a-2.0E0*3*I_ESP_Fx2y_Py_a+2*1*I_ESP_Px_Py;
    abcd[iGrid*180+104] = 4.0E0*I_ESP_Hx3yz_Py_aa-2.0E0*1*I_ESP_Fxyz_Py_a-2.0E0*2*I_ESP_Fxyz_Py_a;
    abcd[iGrid*180+105] = 4.0E0*I_ESP_Hx2y2z_Py_aa-2.0E0*1*I_ESP_Fx2z_Py_a;
    abcd[iGrid*180+106] = 4.0E0*I_ESP_H5y_Py_aa-2.0E0*3*I_ESP_F3y_Py_a-2.0E0*4*I_ESP_F3y_Py_a+3*2*I_ESP_Py_Py;
    abcd[iGrid*180+107] = 4.0E0*I_ESP_H4yz_Py_aa-2.0E0*2*I_ESP_F2yz_Py_a-2.0E0*3*I_ESP_F2yz_Py_a+2*1*I_ESP_Pz_Py;
    abcd[iGrid*180+108] = 4.0E0*I_ESP_H3y2z_Py_aa-2.0E0*1*I_ESP_Fy2z_Py_a-2.0E0*2*I_ESP_Fy2z_Py_a;
    abcd[iGrid*180+109] = 4.0E0*I_ESP_H2y3z_Py_aa-2.0E0*1*I_ESP_F3z_Py_a;
    abcd[iGrid*180+110] = 4.0E0*I_ESP_H3x2y_Pz_aa-2.0E0*1*I_ESP_F3x_Pz_a;
    abcd[iGrid*180+111] = 4.0E0*I_ESP_H2x3y_Pz_aa-2.0E0*1*I_ESP_F2xy_Pz_a-2.0E0*2*I_ESP_F2xy_Pz_a;
    abcd[iGrid*180+112] = 4.0E0*I_ESP_H2x2yz_Pz_aa-2.0E0*1*I_ESP_F2xz_Pz_a;
    abcd[iGrid*180+113] = 4.0E0*I_ESP_Hx4y_Pz_aa-2.0E0*2*I_ESP_Fx2y_Pz_a-2.0E0*3*I_ESP_Fx2y_Pz_a+2*1*I_ESP_Px_Pz;
    abcd[iGrid*180+114] = 4.0E0*I_ESP_Hx3yz_Pz_aa-2.0E0*1*I_ESP_Fxyz_Pz_a-2.0E0*2*I_ESP_Fxyz_Pz_a;
    abcd[iGrid*180+115] = 4.0E0*I_ESP_Hx2y2z_Pz_aa-2.0E0*1*I_ESP_Fx2z_Pz_a;
    abcd[iGrid*180+116] = 4.0E0*I_ESP_H5y_Pz_aa-2.0E0*3*I_ESP_F3y_Pz_a-2.0E0*4*I_ESP_F3y_Pz_a+3*2*I_ESP_Py_Pz;
    abcd[iGrid*180+117] = 4.0E0*I_ESP_H4yz_Pz_aa-2.0E0*2*I_ESP_F2yz_Pz_a-2.0E0*3*I_ESP_F2yz_Pz_a+2*1*I_ESP_Pz_Pz;
    abcd[iGrid*180+118] = 4.0E0*I_ESP_H3y2z_Pz_aa-2.0E0*1*I_ESP_Fy2z_Pz_a-2.0E0*2*I_ESP_Fy2z_Pz_a;
    abcd[iGrid*180+119] = 4.0E0*I_ESP_H2y3z_Pz_aa-2.0E0*1*I_ESP_F3z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_P_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_P_aa
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_P_P
     ************************************************************/
    abcd[iGrid*180+120] = 4.0E0*I_ESP_H3xyz_Px_aa;
    abcd[iGrid*180+121] = 4.0E0*I_ESP_H2x2yz_Px_aa-2.0E0*1*I_ESP_F2xz_Px_a;
    abcd[iGrid*180+122] = 4.0E0*I_ESP_H2xy2z_Px_aa-2.0E0*1*I_ESP_F2xy_Px_a;
    abcd[iGrid*180+123] = 4.0E0*I_ESP_Hx3yz_Px_aa-2.0E0*2*I_ESP_Fxyz_Px_a;
    abcd[iGrid*180+124] = 4.0E0*I_ESP_Hx2y2z_Px_aa-2.0E0*1*I_ESP_Fx2y_Px_a-2.0E0*1*I_ESP_Fx2z_Px_a+1*I_ESP_Px_Px;
    abcd[iGrid*180+125] = 4.0E0*I_ESP_Hxy3z_Px_aa-2.0E0*2*I_ESP_Fxyz_Px_a;
    abcd[iGrid*180+126] = 4.0E0*I_ESP_H4yz_Px_aa-2.0E0*3*I_ESP_F2yz_Px_a;
    abcd[iGrid*180+127] = 4.0E0*I_ESP_H3y2z_Px_aa-2.0E0*1*I_ESP_F3y_Px_a-2.0E0*2*I_ESP_Fy2z_Px_a+2*1*I_ESP_Py_Px;
    abcd[iGrid*180+128] = 4.0E0*I_ESP_H2y3z_Px_aa-2.0E0*2*I_ESP_F2yz_Px_a-2.0E0*1*I_ESP_F3z_Px_a+2*I_ESP_Pz_Px;
    abcd[iGrid*180+129] = 4.0E0*I_ESP_Hy4z_Px_aa-2.0E0*3*I_ESP_Fy2z_Px_a;
    abcd[iGrid*180+130] = 4.0E0*I_ESP_H3xyz_Py_aa;
    abcd[iGrid*180+131] = 4.0E0*I_ESP_H2x2yz_Py_aa-2.0E0*1*I_ESP_F2xz_Py_a;
    abcd[iGrid*180+132] = 4.0E0*I_ESP_H2xy2z_Py_aa-2.0E0*1*I_ESP_F2xy_Py_a;
    abcd[iGrid*180+133] = 4.0E0*I_ESP_Hx3yz_Py_aa-2.0E0*2*I_ESP_Fxyz_Py_a;
    abcd[iGrid*180+134] = 4.0E0*I_ESP_Hx2y2z_Py_aa-2.0E0*1*I_ESP_Fx2y_Py_a-2.0E0*1*I_ESP_Fx2z_Py_a+1*I_ESP_Px_Py;
    abcd[iGrid*180+135] = 4.0E0*I_ESP_Hxy3z_Py_aa-2.0E0*2*I_ESP_Fxyz_Py_a;
    abcd[iGrid*180+136] = 4.0E0*I_ESP_H4yz_Py_aa-2.0E0*3*I_ESP_F2yz_Py_a;
    abcd[iGrid*180+137] = 4.0E0*I_ESP_H3y2z_Py_aa-2.0E0*1*I_ESP_F3y_Py_a-2.0E0*2*I_ESP_Fy2z_Py_a+2*1*I_ESP_Py_Py;
    abcd[iGrid*180+138] = 4.0E0*I_ESP_H2y3z_Py_aa-2.0E0*2*I_ESP_F2yz_Py_a-2.0E0*1*I_ESP_F3z_Py_a+2*I_ESP_Pz_Py;
    abcd[iGrid*180+139] = 4.0E0*I_ESP_Hy4z_Py_aa-2.0E0*3*I_ESP_Fy2z_Py_a;
    abcd[iGrid*180+140] = 4.0E0*I_ESP_H3xyz_Pz_aa;
    abcd[iGrid*180+141] = 4.0E0*I_ESP_H2x2yz_Pz_aa-2.0E0*1*I_ESP_F2xz_Pz_a;
    abcd[iGrid*180+142] = 4.0E0*I_ESP_H2xy2z_Pz_aa-2.0E0*1*I_ESP_F2xy_Pz_a;
    abcd[iGrid*180+143] = 4.0E0*I_ESP_Hx3yz_Pz_aa-2.0E0*2*I_ESP_Fxyz_Pz_a;
    abcd[iGrid*180+144] = 4.0E0*I_ESP_Hx2y2z_Pz_aa-2.0E0*1*I_ESP_Fx2y_Pz_a-2.0E0*1*I_ESP_Fx2z_Pz_a+1*I_ESP_Px_Pz;
    abcd[iGrid*180+145] = 4.0E0*I_ESP_Hxy3z_Pz_aa-2.0E0*2*I_ESP_Fxyz_Pz_a;
    abcd[iGrid*180+146] = 4.0E0*I_ESP_H4yz_Pz_aa-2.0E0*3*I_ESP_F2yz_Pz_a;
    abcd[iGrid*180+147] = 4.0E0*I_ESP_H3y2z_Pz_aa-2.0E0*1*I_ESP_F3y_Pz_a-2.0E0*2*I_ESP_Fy2z_Pz_a+2*1*I_ESP_Py_Pz;
    abcd[iGrid*180+148] = 4.0E0*I_ESP_H2y3z_Pz_aa-2.0E0*2*I_ESP_F2yz_Pz_a-2.0E0*1*I_ESP_F3z_Pz_a+2*I_ESP_Pz_Pz;
    abcd[iGrid*180+149] = 4.0E0*I_ESP_Hy4z_Pz_aa-2.0E0*3*I_ESP_Fy2z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_P_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_P_aa
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_F_P_a
     * RHS shell quartet name: SQ_ESP_P_P
     ************************************************************/
    abcd[iGrid*180+150] = 4.0E0*I_ESP_H3x2z_Px_aa-2.0E0*1*I_ESP_F3x_Px_a;
    abcd[iGrid*180+151] = 4.0E0*I_ESP_H2xy2z_Px_aa-2.0E0*1*I_ESP_F2xy_Px_a;
    abcd[iGrid*180+152] = 4.0E0*I_ESP_H2x3z_Px_aa-2.0E0*1*I_ESP_F2xz_Px_a-2.0E0*2*I_ESP_F2xz_Px_a;
    abcd[iGrid*180+153] = 4.0E0*I_ESP_Hx2y2z_Px_aa-2.0E0*1*I_ESP_Fx2y_Px_a;
    abcd[iGrid*180+154] = 4.0E0*I_ESP_Hxy3z_Px_aa-2.0E0*1*I_ESP_Fxyz_Px_a-2.0E0*2*I_ESP_Fxyz_Px_a;
    abcd[iGrid*180+155] = 4.0E0*I_ESP_Hx4z_Px_aa-2.0E0*2*I_ESP_Fx2z_Px_a-2.0E0*3*I_ESP_Fx2z_Px_a+2*1*I_ESP_Px_Px;
    abcd[iGrid*180+156] = 4.0E0*I_ESP_H3y2z_Px_aa-2.0E0*1*I_ESP_F3y_Px_a;
    abcd[iGrid*180+157] = 4.0E0*I_ESP_H2y3z_Px_aa-2.0E0*1*I_ESP_F2yz_Px_a-2.0E0*2*I_ESP_F2yz_Px_a;
    abcd[iGrid*180+158] = 4.0E0*I_ESP_Hy4z_Px_aa-2.0E0*2*I_ESP_Fy2z_Px_a-2.0E0*3*I_ESP_Fy2z_Px_a+2*1*I_ESP_Py_Px;
    abcd[iGrid*180+159] = 4.0E0*I_ESP_H5z_Px_aa-2.0E0*3*I_ESP_F3z_Px_a-2.0E0*4*I_ESP_F3z_Px_a+3*2*I_ESP_Pz_Px;
    abcd[iGrid*180+160] = 4.0E0*I_ESP_H3x2z_Py_aa-2.0E0*1*I_ESP_F3x_Py_a;
    abcd[iGrid*180+161] = 4.0E0*I_ESP_H2xy2z_Py_aa-2.0E0*1*I_ESP_F2xy_Py_a;
    abcd[iGrid*180+162] = 4.0E0*I_ESP_H2x3z_Py_aa-2.0E0*1*I_ESP_F2xz_Py_a-2.0E0*2*I_ESP_F2xz_Py_a;
    abcd[iGrid*180+163] = 4.0E0*I_ESP_Hx2y2z_Py_aa-2.0E0*1*I_ESP_Fx2y_Py_a;
    abcd[iGrid*180+164] = 4.0E0*I_ESP_Hxy3z_Py_aa-2.0E0*1*I_ESP_Fxyz_Py_a-2.0E0*2*I_ESP_Fxyz_Py_a;
    abcd[iGrid*180+165] = 4.0E0*I_ESP_Hx4z_Py_aa-2.0E0*2*I_ESP_Fx2z_Py_a-2.0E0*3*I_ESP_Fx2z_Py_a+2*1*I_ESP_Px_Py;
    abcd[iGrid*180+166] = 4.0E0*I_ESP_H3y2z_Py_aa-2.0E0*1*I_ESP_F3y_Py_a;
    abcd[iGrid*180+167] = 4.0E0*I_ESP_H2y3z_Py_aa-2.0E0*1*I_ESP_F2yz_Py_a-2.0E0*2*I_ESP_F2yz_Py_a;
    abcd[iGrid*180+168] = 4.0E0*I_ESP_Hy4z_Py_aa-2.0E0*2*I_ESP_Fy2z_Py_a-2.0E0*3*I_ESP_Fy2z_Py_a+2*1*I_ESP_Py_Py;
    abcd[iGrid*180+169] = 4.0E0*I_ESP_H5z_Py_aa-2.0E0*3*I_ESP_F3z_Py_a-2.0E0*4*I_ESP_F3z_Py_a+3*2*I_ESP_Pz_Py;
    abcd[iGrid*180+170] = 4.0E0*I_ESP_H3x2z_Pz_aa-2.0E0*1*I_ESP_F3x_Pz_a;
    abcd[iGrid*180+171] = 4.0E0*I_ESP_H2xy2z_Pz_aa-2.0E0*1*I_ESP_F2xy_Pz_a;
    abcd[iGrid*180+172] = 4.0E0*I_ESP_H2x3z_Pz_aa-2.0E0*1*I_ESP_F2xz_Pz_a-2.0E0*2*I_ESP_F2xz_Pz_a;
    abcd[iGrid*180+173] = 4.0E0*I_ESP_Hx2y2z_Pz_aa-2.0E0*1*I_ESP_Fx2y_Pz_a;
    abcd[iGrid*180+174] = 4.0E0*I_ESP_Hxy3z_Pz_aa-2.0E0*1*I_ESP_Fxyz_Pz_a-2.0E0*2*I_ESP_Fxyz_Pz_a;
    abcd[iGrid*180+175] = 4.0E0*I_ESP_Hx4z_Pz_aa-2.0E0*2*I_ESP_Fx2z_Pz_a-2.0E0*3*I_ESP_Fx2z_Pz_a+2*1*I_ESP_Px_Pz;
    abcd[iGrid*180+176] = 4.0E0*I_ESP_H3y2z_Pz_aa-2.0E0*1*I_ESP_F3y_Pz_a;
    abcd[iGrid*180+177] = 4.0E0*I_ESP_H2y3z_Pz_aa-2.0E0*1*I_ESP_F2yz_Pz_a-2.0E0*2*I_ESP_F2yz_Pz_a;
    abcd[iGrid*180+178] = 4.0E0*I_ESP_Hy4z_Pz_aa-2.0E0*2*I_ESP_Fy2z_Pz_a-2.0E0*3*I_ESP_Fy2z_Pz_a+2*1*I_ESP_Py_Pz;
    abcd[iGrid*180+179] = 4.0E0*I_ESP_H5z_Pz_aa-2.0E0*3*I_ESP_F3z_Pz_a-2.0E0*4*I_ESP_F3z_Pz_a+3*2*I_ESP_Pz_Pz;
  }
}
