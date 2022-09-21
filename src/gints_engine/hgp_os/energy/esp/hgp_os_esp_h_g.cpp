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

void hgp_os_esp_h_g(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_M9x_S = 0.0E0;
    Double I_ESP_M8xy_S = 0.0E0;
    Double I_ESP_M8xz_S = 0.0E0;
    Double I_ESP_M7x2y_S = 0.0E0;
    Double I_ESP_M7xyz_S = 0.0E0;
    Double I_ESP_M7x2z_S = 0.0E0;
    Double I_ESP_M6x3y_S = 0.0E0;
    Double I_ESP_M6x2yz_S = 0.0E0;
    Double I_ESP_M6xy2z_S = 0.0E0;
    Double I_ESP_M6x3z_S = 0.0E0;
    Double I_ESP_M5x4y_S = 0.0E0;
    Double I_ESP_M5x3yz_S = 0.0E0;
    Double I_ESP_M5x2y2z_S = 0.0E0;
    Double I_ESP_M5xy3z_S = 0.0E0;
    Double I_ESP_M5x4z_S = 0.0E0;
    Double I_ESP_M4x5y_S = 0.0E0;
    Double I_ESP_M4x4yz_S = 0.0E0;
    Double I_ESP_M4x3y2z_S = 0.0E0;
    Double I_ESP_M4x2y3z_S = 0.0E0;
    Double I_ESP_M4xy4z_S = 0.0E0;
    Double I_ESP_M4x5z_S = 0.0E0;
    Double I_ESP_M3x6y_S = 0.0E0;
    Double I_ESP_M3x5yz_S = 0.0E0;
    Double I_ESP_M3x4y2z_S = 0.0E0;
    Double I_ESP_M3x3y3z_S = 0.0E0;
    Double I_ESP_M3x2y4z_S = 0.0E0;
    Double I_ESP_M3xy5z_S = 0.0E0;
    Double I_ESP_M3x6z_S = 0.0E0;
    Double I_ESP_M2x7y_S = 0.0E0;
    Double I_ESP_M2x6yz_S = 0.0E0;
    Double I_ESP_M2x5y2z_S = 0.0E0;
    Double I_ESP_M2x4y3z_S = 0.0E0;
    Double I_ESP_M2x3y4z_S = 0.0E0;
    Double I_ESP_M2x2y5z_S = 0.0E0;
    Double I_ESP_M2xy6z_S = 0.0E0;
    Double I_ESP_M2x7z_S = 0.0E0;
    Double I_ESP_Mx8y_S = 0.0E0;
    Double I_ESP_Mx7yz_S = 0.0E0;
    Double I_ESP_Mx6y2z_S = 0.0E0;
    Double I_ESP_Mx5y3z_S = 0.0E0;
    Double I_ESP_Mx4y4z_S = 0.0E0;
    Double I_ESP_Mx3y5z_S = 0.0E0;
    Double I_ESP_Mx2y6z_S = 0.0E0;
    Double I_ESP_Mxy7z_S = 0.0E0;
    Double I_ESP_Mx8z_S = 0.0E0;
    Double I_ESP_M9y_S = 0.0E0;
    Double I_ESP_M8yz_S = 0.0E0;
    Double I_ESP_M7y2z_S = 0.0E0;
    Double I_ESP_M6y3z_S = 0.0E0;
    Double I_ESP_M5y4z_S = 0.0E0;
    Double I_ESP_M4y5z_S = 0.0E0;
    Double I_ESP_M3y6z_S = 0.0E0;
    Double I_ESP_M2y7z_S = 0.0E0;
    Double I_ESP_My8z_S = 0.0E0;
    Double I_ESP_M9z_S = 0.0E0;
    Double I_ESP_L8x_S = 0.0E0;
    Double I_ESP_L7xy_S = 0.0E0;
    Double I_ESP_L7xz_S = 0.0E0;
    Double I_ESP_L6x2y_S = 0.0E0;
    Double I_ESP_L6xyz_S = 0.0E0;
    Double I_ESP_L6x2z_S = 0.0E0;
    Double I_ESP_L5x3y_S = 0.0E0;
    Double I_ESP_L5x2yz_S = 0.0E0;
    Double I_ESP_L5xy2z_S = 0.0E0;
    Double I_ESP_L5x3z_S = 0.0E0;
    Double I_ESP_L4x4y_S = 0.0E0;
    Double I_ESP_L4x3yz_S = 0.0E0;
    Double I_ESP_L4x2y2z_S = 0.0E0;
    Double I_ESP_L4xy3z_S = 0.0E0;
    Double I_ESP_L4x4z_S = 0.0E0;
    Double I_ESP_L3x5y_S = 0.0E0;
    Double I_ESP_L3x4yz_S = 0.0E0;
    Double I_ESP_L3x3y2z_S = 0.0E0;
    Double I_ESP_L3x2y3z_S = 0.0E0;
    Double I_ESP_L3xy4z_S = 0.0E0;
    Double I_ESP_L3x5z_S = 0.0E0;
    Double I_ESP_L2x6y_S = 0.0E0;
    Double I_ESP_L2x5yz_S = 0.0E0;
    Double I_ESP_L2x4y2z_S = 0.0E0;
    Double I_ESP_L2x3y3z_S = 0.0E0;
    Double I_ESP_L2x2y4z_S = 0.0E0;
    Double I_ESP_L2xy5z_S = 0.0E0;
    Double I_ESP_L2x6z_S = 0.0E0;
    Double I_ESP_Lx7y_S = 0.0E0;
    Double I_ESP_Lx6yz_S = 0.0E0;
    Double I_ESP_Lx5y2z_S = 0.0E0;
    Double I_ESP_Lx4y3z_S = 0.0E0;
    Double I_ESP_Lx3y4z_S = 0.0E0;
    Double I_ESP_Lx2y5z_S = 0.0E0;
    Double I_ESP_Lxy6z_S = 0.0E0;
    Double I_ESP_Lx7z_S = 0.0E0;
    Double I_ESP_L8y_S = 0.0E0;
    Double I_ESP_L7yz_S = 0.0E0;
    Double I_ESP_L6y2z_S = 0.0E0;
    Double I_ESP_L5y3z_S = 0.0E0;
    Double I_ESP_L4y4z_S = 0.0E0;
    Double I_ESP_L3y5z_S = 0.0E0;
    Double I_ESP_L2y6z_S = 0.0E0;
    Double I_ESP_Ly7z_S = 0.0E0;
    Double I_ESP_L8z_S = 0.0E0;
    Double I_ESP_K7x_S = 0.0E0;
    Double I_ESP_K6xy_S = 0.0E0;
    Double I_ESP_K6xz_S = 0.0E0;
    Double I_ESP_K5x2y_S = 0.0E0;
    Double I_ESP_K5xyz_S = 0.0E0;
    Double I_ESP_K5x2z_S = 0.0E0;
    Double I_ESP_K4x3y_S = 0.0E0;
    Double I_ESP_K4x2yz_S = 0.0E0;
    Double I_ESP_K4xy2z_S = 0.0E0;
    Double I_ESP_K4x3z_S = 0.0E0;
    Double I_ESP_K3x4y_S = 0.0E0;
    Double I_ESP_K3x3yz_S = 0.0E0;
    Double I_ESP_K3x2y2z_S = 0.0E0;
    Double I_ESP_K3xy3z_S = 0.0E0;
    Double I_ESP_K3x4z_S = 0.0E0;
    Double I_ESP_K2x5y_S = 0.0E0;
    Double I_ESP_K2x4yz_S = 0.0E0;
    Double I_ESP_K2x3y2z_S = 0.0E0;
    Double I_ESP_K2x2y3z_S = 0.0E0;
    Double I_ESP_K2xy4z_S = 0.0E0;
    Double I_ESP_K2x5z_S = 0.0E0;
    Double I_ESP_Kx6y_S = 0.0E0;
    Double I_ESP_Kx5yz_S = 0.0E0;
    Double I_ESP_Kx4y2z_S = 0.0E0;
    Double I_ESP_Kx3y3z_S = 0.0E0;
    Double I_ESP_Kx2y4z_S = 0.0E0;
    Double I_ESP_Kxy5z_S = 0.0E0;
    Double I_ESP_Kx6z_S = 0.0E0;
    Double I_ESP_K7y_S = 0.0E0;
    Double I_ESP_K6yz_S = 0.0E0;
    Double I_ESP_K5y2z_S = 0.0E0;
    Double I_ESP_K4y3z_S = 0.0E0;
    Double I_ESP_K3y4z_S = 0.0E0;
    Double I_ESP_K2y5z_S = 0.0E0;
    Double I_ESP_Ky6z_S = 0.0E0;
    Double I_ESP_K7z_S = 0.0E0;
    Double I_ESP_I6x_S = 0.0E0;
    Double I_ESP_I5xy_S = 0.0E0;
    Double I_ESP_I5xz_S = 0.0E0;
    Double I_ESP_I4x2y_S = 0.0E0;
    Double I_ESP_I4xyz_S = 0.0E0;
    Double I_ESP_I4x2z_S = 0.0E0;
    Double I_ESP_I3x3y_S = 0.0E0;
    Double I_ESP_I3x2yz_S = 0.0E0;
    Double I_ESP_I3xy2z_S = 0.0E0;
    Double I_ESP_I3x3z_S = 0.0E0;
    Double I_ESP_I2x4y_S = 0.0E0;
    Double I_ESP_I2x3yz_S = 0.0E0;
    Double I_ESP_I2x2y2z_S = 0.0E0;
    Double I_ESP_I2xy3z_S = 0.0E0;
    Double I_ESP_I2x4z_S = 0.0E0;
    Double I_ESP_Ix5y_S = 0.0E0;
    Double I_ESP_Ix4yz_S = 0.0E0;
    Double I_ESP_Ix3y2z_S = 0.0E0;
    Double I_ESP_Ix2y3z_S = 0.0E0;
    Double I_ESP_Ixy4z_S = 0.0E0;
    Double I_ESP_Ix5z_S = 0.0E0;
    Double I_ESP_I6y_S = 0.0E0;
    Double I_ESP_I5yz_S = 0.0E0;
    Double I_ESP_I4y2z_S = 0.0E0;
    Double I_ESP_I3y3z_S = 0.0E0;
    Double I_ESP_I2y4z_S = 0.0E0;
    Double I_ESP_Iy5z_S = 0.0E0;
    Double I_ESP_I6z_S = 0.0E0;
    Double I_ESP_H5x_S = 0.0E0;
    Double I_ESP_H4xy_S = 0.0E0;
    Double I_ESP_H4xz_S = 0.0E0;
    Double I_ESP_H3x2y_S = 0.0E0;
    Double I_ESP_H3xyz_S = 0.0E0;
    Double I_ESP_H3x2z_S = 0.0E0;
    Double I_ESP_H2x3y_S = 0.0E0;
    Double I_ESP_H2x2yz_S = 0.0E0;
    Double I_ESP_H2xy2z_S = 0.0E0;
    Double I_ESP_H2x3z_S = 0.0E0;
    Double I_ESP_Hx4y_S = 0.0E0;
    Double I_ESP_Hx3yz_S = 0.0E0;
    Double I_ESP_Hx2y2z_S = 0.0E0;
    Double I_ESP_Hxy3z_S = 0.0E0;
    Double I_ESP_Hx4z_S = 0.0E0;
    Double I_ESP_H5y_S = 0.0E0;
    Double I_ESP_H4yz_S = 0.0E0;
    Double I_ESP_H3y2z_S = 0.0E0;
    Double I_ESP_H2y3z_S = 0.0E0;
    Double I_ESP_Hy4z_S = 0.0E0;
    Double I_ESP_H5z_S = 0.0E0;

    for(UInt ip2=0; ip2<inp2; ip2++) {
      Double ic2   = icoe[ip2];
      Double onedz = iexp[ip2];
      Double rho   = 1.0E0/onedz;
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
      Double I_ESP_S_S_M7_vrr  = 0.0E0;
      Double I_ESP_S_S_M8_vrr  = 0.0E0;
      Double I_ESP_S_S_M9_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ESP_S_S_vrr_d  = 0.0E0;
      double I_ESP_S_S_M1_vrr_d  = 0.0E0;
      double I_ESP_S_S_M2_vrr_d  = 0.0E0;
      double I_ESP_S_S_M3_vrr_d  = 0.0E0;
      double I_ESP_S_S_M4_vrr_d  = 0.0E0;
      double I_ESP_S_S_M5_vrr_d  = 0.0E0;
      double I_ESP_S_S_M6_vrr_d  = 0.0E0;
      double I_ESP_S_S_M7_vrr_d  = 0.0E0;
      double I_ESP_S_S_M8_vrr_d  = 0.0E0;
      double I_ESP_S_S_M9_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER53;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER51*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER49*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER47*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER45*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER43*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER41*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER39*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER37*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER35*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER33*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER31*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER29*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER27*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER25*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER23*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER21*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = ONEOVER19*I_ESP_S_S_M9_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M9_vrr  = f*I_ESP_S_S_M9_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ESP_S_S_M8_vrr  = ONEOVER17*(u2*I_ESP_S_S_M9_vrr+f);
        I_ESP_S_S_M7_vrr  = ONEOVER15*(u2*I_ESP_S_S_M8_vrr+f);
        I_ESP_S_S_M6_vrr  = ONEOVER13*(u2*I_ESP_S_S_M7_vrr+f);
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
        I_ESP_S_S_M7_vrr_d = oneO2u*(13.0E0*I_ESP_S_S_M6_vrr_d-f);
        I_ESP_S_S_M8_vrr_d = oneO2u*(15.0E0*I_ESP_S_S_M7_vrr_d-f);
        I_ESP_S_S_M9_vrr_d = oneO2u*(17.0E0*I_ESP_S_S_M8_vrr_d-f);

        // write the double result back to the float var
        I_ESP_S_S_vrr = static_cast<Double>(I_ESP_S_S_vrr_d);
        I_ESP_S_S_M1_vrr = static_cast<Double>(I_ESP_S_S_M1_vrr_d);
        I_ESP_S_S_M2_vrr = static_cast<Double>(I_ESP_S_S_M2_vrr_d);
        I_ESP_S_S_M3_vrr = static_cast<Double>(I_ESP_S_S_M3_vrr_d);
        I_ESP_S_S_M4_vrr = static_cast<Double>(I_ESP_S_S_M4_vrr_d);
        I_ESP_S_S_M5_vrr = static_cast<Double>(I_ESP_S_S_M5_vrr_d);
        I_ESP_S_S_M6_vrr = static_cast<Double>(I_ESP_S_S_M6_vrr_d);
        I_ESP_S_S_M7_vrr = static_cast<Double>(I_ESP_S_S_M7_vrr_d);
        I_ESP_S_S_M8_vrr = static_cast<Double>(I_ESP_S_S_M8_vrr_d);
        I_ESP_S_S_M9_vrr = static_cast<Double>(I_ESP_S_S_M9_vrr_d);

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
        I_ESP_S_S_M7_vrr = oneO2u*(13.0E0*I_ESP_S_S_M6_vrr-f);
        I_ESP_S_S_M8_vrr = oneO2u*(15.0E0*I_ESP_S_S_M7_vrr-f);
        I_ESP_S_S_M9_vrr = oneO2u*(17.0E0*I_ESP_S_S_M8_vrr-f);

#endif

      }


        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M8
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M8
         * RHS shell quartet name: SQ_ESP_S_S_M9
         ************************************************************/
        Double I_ESP_Px_S_M8_vrr = PAX*I_ESP_S_S_M8_vrr-PRX*I_ESP_S_S_M9_vrr;
        Double I_ESP_Py_S_M8_vrr = PAY*I_ESP_S_S_M8_vrr-PRY*I_ESP_S_S_M9_vrr;
        Double I_ESP_Pz_S_M8_vrr = PAZ*I_ESP_S_S_M8_vrr-PRZ*I_ESP_S_S_M9_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M7
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M7
         * RHS shell quartet name: SQ_ESP_S_S_M8
         ************************************************************/
        Double I_ESP_Px_S_M7_vrr = PAX*I_ESP_S_S_M7_vrr-PRX*I_ESP_S_S_M8_vrr;
        Double I_ESP_Py_S_M7_vrr = PAY*I_ESP_S_S_M7_vrr-PRY*I_ESP_S_S_M8_vrr;
        Double I_ESP_Pz_S_M7_vrr = PAZ*I_ESP_S_S_M7_vrr-PRZ*I_ESP_S_S_M8_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M7
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M7
         * RHS shell quartet name: SQ_ESP_P_S_M8
         * RHS shell quartet name: SQ_ESP_S_S_M7
         * RHS shell quartet name: SQ_ESP_S_S_M8
         ************************************************************/
        Double I_ESP_D2x_S_M7_vrr = PAX*I_ESP_Px_S_M7_vrr-PRX*I_ESP_Px_S_M8_vrr+oned2z*I_ESP_S_S_M7_vrr-oned2z*I_ESP_S_S_M8_vrr;
        Double I_ESP_D2y_S_M7_vrr = PAY*I_ESP_Py_S_M7_vrr-PRY*I_ESP_Py_S_M8_vrr+oned2z*I_ESP_S_S_M7_vrr-oned2z*I_ESP_S_S_M8_vrr;
        Double I_ESP_D2z_S_M7_vrr = PAZ*I_ESP_Pz_S_M7_vrr-PRZ*I_ESP_Pz_S_M8_vrr+oned2z*I_ESP_S_S_M7_vrr-oned2z*I_ESP_S_S_M8_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M6
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M6
         * RHS shell quartet name: SQ_ESP_S_S_M7
         ************************************************************/
        Double I_ESP_Px_S_M6_vrr = PAX*I_ESP_S_S_M6_vrr-PRX*I_ESP_S_S_M7_vrr;
        Double I_ESP_Py_S_M6_vrr = PAY*I_ESP_S_S_M6_vrr-PRY*I_ESP_S_S_M7_vrr;
        Double I_ESP_Pz_S_M6_vrr = PAZ*I_ESP_S_S_M6_vrr-PRZ*I_ESP_S_S_M7_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M6
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M6
         * RHS shell quartet name: SQ_ESP_P_S_M7
         * RHS shell quartet name: SQ_ESP_S_S_M6
         * RHS shell quartet name: SQ_ESP_S_S_M7
         ************************************************************/
        Double I_ESP_D2x_S_M6_vrr = PAX*I_ESP_Px_S_M6_vrr-PRX*I_ESP_Px_S_M7_vrr+oned2z*I_ESP_S_S_M6_vrr-oned2z*I_ESP_S_S_M7_vrr;
        Double I_ESP_D2y_S_M6_vrr = PAY*I_ESP_Py_S_M6_vrr-PRY*I_ESP_Py_S_M7_vrr+oned2z*I_ESP_S_S_M6_vrr-oned2z*I_ESP_S_S_M7_vrr;
        Double I_ESP_D2z_S_M6_vrr = PAZ*I_ESP_Pz_S_M6_vrr-PRZ*I_ESP_Pz_S_M7_vrr+oned2z*I_ESP_S_S_M6_vrr-oned2z*I_ESP_S_S_M7_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M6
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M6
         * RHS shell quartet name: SQ_ESP_D_S_M7
         * RHS shell quartet name: SQ_ESP_P_S_M6
         * RHS shell quartet name: SQ_ESP_P_S_M7
         ************************************************************/
        Double I_ESP_F3x_S_M6_vrr = PAX*I_ESP_D2x_S_M6_vrr-PRX*I_ESP_D2x_S_M7_vrr+2*oned2z*I_ESP_Px_S_M6_vrr-2*oned2z*I_ESP_Px_S_M7_vrr;
        Double I_ESP_F3y_S_M6_vrr = PAY*I_ESP_D2y_S_M6_vrr-PRY*I_ESP_D2y_S_M7_vrr+2*oned2z*I_ESP_Py_S_M6_vrr-2*oned2z*I_ESP_Py_S_M7_vrr;
        Double I_ESP_F3z_S_M6_vrr = PAZ*I_ESP_D2z_S_M6_vrr-PRZ*I_ESP_D2z_S_M7_vrr+2*oned2z*I_ESP_Pz_S_M6_vrr-2*oned2z*I_ESP_Pz_S_M7_vrr;

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
         * shell quartet name: SQ_ESP_D_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M5
         * RHS shell quartet name: SQ_ESP_P_S_M6
         * RHS shell quartet name: SQ_ESP_S_S_M5
         * RHS shell quartet name: SQ_ESP_S_S_M6
         ************************************************************/
        Double I_ESP_D2x_S_M5_vrr = PAX*I_ESP_Px_S_M5_vrr-PRX*I_ESP_Px_S_M6_vrr+oned2z*I_ESP_S_S_M5_vrr-oned2z*I_ESP_S_S_M6_vrr;
        Double I_ESP_D2y_S_M5_vrr = PAY*I_ESP_Py_S_M5_vrr-PRY*I_ESP_Py_S_M6_vrr+oned2z*I_ESP_S_S_M5_vrr-oned2z*I_ESP_S_S_M6_vrr;
        Double I_ESP_D2z_S_M5_vrr = PAZ*I_ESP_Pz_S_M5_vrr-PRZ*I_ESP_Pz_S_M6_vrr+oned2z*I_ESP_S_S_M5_vrr-oned2z*I_ESP_S_S_M6_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M5
         * RHS shell quartet name: SQ_ESP_D_S_M6
         * RHS shell quartet name: SQ_ESP_P_S_M5
         * RHS shell quartet name: SQ_ESP_P_S_M6
         ************************************************************/
        Double I_ESP_F3x_S_M5_vrr = PAX*I_ESP_D2x_S_M5_vrr-PRX*I_ESP_D2x_S_M6_vrr+2*oned2z*I_ESP_Px_S_M5_vrr-2*oned2z*I_ESP_Px_S_M6_vrr;
        Double I_ESP_F3y_S_M5_vrr = PAY*I_ESP_D2y_S_M5_vrr-PRY*I_ESP_D2y_S_M6_vrr+2*oned2z*I_ESP_Py_S_M5_vrr-2*oned2z*I_ESP_Py_S_M6_vrr;
        Double I_ESP_F3z_S_M5_vrr = PAZ*I_ESP_D2z_S_M5_vrr-PRZ*I_ESP_D2z_S_M6_vrr+2*oned2z*I_ESP_Pz_S_M5_vrr-2*oned2z*I_ESP_Pz_S_M6_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 11 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M5
         * RHS shell quartet name: SQ_ESP_F_S_M6
         * RHS shell quartet name: SQ_ESP_D_S_M5
         * RHS shell quartet name: SQ_ESP_D_S_M6
         ************************************************************/
        Double I_ESP_G4x_S_M5_vrr = PAX*I_ESP_F3x_S_M5_vrr-PRX*I_ESP_F3x_S_M6_vrr+3*oned2z*I_ESP_D2x_S_M5_vrr-3*oned2z*I_ESP_D2x_S_M6_vrr;
        Double I_ESP_G3xy_S_M5_vrr = PAY*I_ESP_F3x_S_M5_vrr-PRY*I_ESP_F3x_S_M6_vrr;
        Double I_ESP_G4y_S_M5_vrr = PAY*I_ESP_F3y_S_M5_vrr-PRY*I_ESP_F3y_S_M6_vrr+3*oned2z*I_ESP_D2y_S_M5_vrr-3*oned2z*I_ESP_D2y_S_M6_vrr;
        Double I_ESP_G4z_S_M5_vrr = PAZ*I_ESP_F3z_S_M5_vrr-PRZ*I_ESP_F3z_S_M6_vrr+3*oned2z*I_ESP_D2z_S_M5_vrr-3*oned2z*I_ESP_D2z_S_M6_vrr;

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
         * shell quartet name: SQ_ESP_F_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M4
         * RHS shell quartet name: SQ_ESP_D_S_M5
         * RHS shell quartet name: SQ_ESP_P_S_M4
         * RHS shell quartet name: SQ_ESP_P_S_M5
         ************************************************************/
        Double I_ESP_F3x_S_M4_vrr = PAX*I_ESP_D2x_S_M4_vrr-PRX*I_ESP_D2x_S_M5_vrr+2*oned2z*I_ESP_Px_S_M4_vrr-2*oned2z*I_ESP_Px_S_M5_vrr;
        Double I_ESP_F3y_S_M4_vrr = PAY*I_ESP_D2y_S_M4_vrr-PRY*I_ESP_D2y_S_M5_vrr+2*oned2z*I_ESP_Py_S_M4_vrr-2*oned2z*I_ESP_Py_S_M5_vrr;
        Double I_ESP_F3z_S_M4_vrr = PAZ*I_ESP_D2z_S_M4_vrr-PRZ*I_ESP_D2z_S_M5_vrr+2*oned2z*I_ESP_Pz_S_M4_vrr-2*oned2z*I_ESP_Pz_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 9 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M4
         * RHS shell quartet name: SQ_ESP_F_S_M5
         * RHS shell quartet name: SQ_ESP_D_S_M4
         * RHS shell quartet name: SQ_ESP_D_S_M5
         ************************************************************/
        Double I_ESP_G4x_S_M4_vrr = PAX*I_ESP_F3x_S_M4_vrr-PRX*I_ESP_F3x_S_M5_vrr+3*oned2z*I_ESP_D2x_S_M4_vrr-3*oned2z*I_ESP_D2x_S_M5_vrr;
        Double I_ESP_G3xy_S_M4_vrr = PAY*I_ESP_F3x_S_M4_vrr-PRY*I_ESP_F3x_S_M5_vrr;
        Double I_ESP_G3xz_S_M4_vrr = PAZ*I_ESP_F3x_S_M4_vrr-PRZ*I_ESP_F3x_S_M5_vrr;
        Double I_ESP_G4y_S_M4_vrr = PAY*I_ESP_F3y_S_M4_vrr-PRY*I_ESP_F3y_S_M5_vrr+3*oned2z*I_ESP_D2y_S_M4_vrr-3*oned2z*I_ESP_D2y_S_M5_vrr;
        Double I_ESP_G3yz_S_M4_vrr = PAZ*I_ESP_F3y_S_M4_vrr-PRZ*I_ESP_F3y_S_M5_vrr;
        Double I_ESP_G4z_S_M4_vrr = PAZ*I_ESP_F3z_S_M4_vrr-PRZ*I_ESP_F3z_S_M5_vrr+3*oned2z*I_ESP_D2z_S_M4_vrr-3*oned2z*I_ESP_D2z_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 11 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M4
         * RHS shell quartet name: SQ_ESP_G_S_M5
         * RHS shell quartet name: SQ_ESP_F_S_M4
         * RHS shell quartet name: SQ_ESP_F_S_M5
         ************************************************************/
        Double I_ESP_H5x_S_M4_vrr = PAX*I_ESP_G4x_S_M4_vrr-PRX*I_ESP_G4x_S_M5_vrr+4*oned2z*I_ESP_F3x_S_M4_vrr-4*oned2z*I_ESP_F3x_S_M5_vrr;
        Double I_ESP_H4xy_S_M4_vrr = PAY*I_ESP_G4x_S_M4_vrr-PRY*I_ESP_G4x_S_M5_vrr;
        Double I_ESP_H4xz_S_M4_vrr = PAZ*I_ESP_G4x_S_M4_vrr-PRZ*I_ESP_G4x_S_M5_vrr;
        Double I_ESP_H3x2y_S_M4_vrr = PAY*I_ESP_G3xy_S_M4_vrr-PRY*I_ESP_G3xy_S_M5_vrr+oned2z*I_ESP_F3x_S_M4_vrr-oned2z*I_ESP_F3x_S_M5_vrr;
        Double I_ESP_Hx4y_S_M4_vrr = PAX*I_ESP_G4y_S_M4_vrr-PRX*I_ESP_G4y_S_M5_vrr;
        Double I_ESP_Hx4z_S_M4_vrr = PAX*I_ESP_G4z_S_M4_vrr-PRX*I_ESP_G4z_S_M5_vrr;
        Double I_ESP_H5y_S_M4_vrr = PAY*I_ESP_G4y_S_M4_vrr-PRY*I_ESP_G4y_S_M5_vrr+4*oned2z*I_ESP_F3y_S_M4_vrr-4*oned2z*I_ESP_F3y_S_M5_vrr;
        Double I_ESP_H4yz_S_M4_vrr = PAZ*I_ESP_G4y_S_M4_vrr-PRZ*I_ESP_G4y_S_M5_vrr;
        Double I_ESP_Hy4z_S_M4_vrr = PAY*I_ESP_G4z_S_M4_vrr-PRY*I_ESP_G4z_S_M5_vrr;
        Double I_ESP_H5z_S_M4_vrr = PAZ*I_ESP_G4z_S_M4_vrr-PRZ*I_ESP_G4z_S_M5_vrr+4*oned2z*I_ESP_F3z_S_M4_vrr-4*oned2z*I_ESP_F3z_S_M5_vrr;

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
         * shell quartet name: SQ_ESP_G_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M3
         * RHS shell quartet name: SQ_ESP_F_S_M4
         * RHS shell quartet name: SQ_ESP_D_S_M3
         * RHS shell quartet name: SQ_ESP_D_S_M4
         ************************************************************/
        Double I_ESP_G4x_S_M3_vrr = PAX*I_ESP_F3x_S_M3_vrr-PRX*I_ESP_F3x_S_M4_vrr+3*oned2z*I_ESP_D2x_S_M3_vrr-3*oned2z*I_ESP_D2x_S_M4_vrr;
        Double I_ESP_G3xy_S_M3_vrr = PAY*I_ESP_F3x_S_M3_vrr-PRY*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_G3xz_S_M3_vrr = PAZ*I_ESP_F3x_S_M3_vrr-PRZ*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_Gx3y_S_M3_vrr = PAX*I_ESP_F3y_S_M3_vrr-PRX*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_Gx3z_S_M3_vrr = PAX*I_ESP_F3z_S_M3_vrr-PRX*I_ESP_F3z_S_M4_vrr;
        Double I_ESP_G4y_S_M3_vrr = PAY*I_ESP_F3y_S_M3_vrr-PRY*I_ESP_F3y_S_M4_vrr+3*oned2z*I_ESP_D2y_S_M3_vrr-3*oned2z*I_ESP_D2y_S_M4_vrr;
        Double I_ESP_G3yz_S_M3_vrr = PAZ*I_ESP_F3y_S_M3_vrr-PRZ*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_G4z_S_M3_vrr = PAZ*I_ESP_F3z_S_M3_vrr-PRZ*I_ESP_F3z_S_M4_vrr+3*oned2z*I_ESP_D2z_S_M3_vrr-3*oned2z*I_ESP_D2z_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 9 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M3
         * RHS shell quartet name: SQ_ESP_G_S_M4
         * RHS shell quartet name: SQ_ESP_F_S_M3
         * RHS shell quartet name: SQ_ESP_F_S_M4
         ************************************************************/
        Double I_ESP_H5x_S_M3_vrr = PAX*I_ESP_G4x_S_M3_vrr-PRX*I_ESP_G4x_S_M4_vrr+4*oned2z*I_ESP_F3x_S_M3_vrr-4*oned2z*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_H4xy_S_M3_vrr = PAY*I_ESP_G4x_S_M3_vrr-PRY*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_H4xz_S_M3_vrr = PAZ*I_ESP_G4x_S_M3_vrr-PRZ*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_H3x2y_S_M3_vrr = PAY*I_ESP_G3xy_S_M3_vrr-PRY*I_ESP_G3xy_S_M4_vrr+oned2z*I_ESP_F3x_S_M3_vrr-oned2z*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_H3x2z_S_M3_vrr = PAZ*I_ESP_G3xz_S_M3_vrr-PRZ*I_ESP_G3xz_S_M4_vrr+oned2z*I_ESP_F3x_S_M3_vrr-oned2z*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_Hx4y_S_M3_vrr = PAX*I_ESP_G4y_S_M3_vrr-PRX*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_Hx4z_S_M3_vrr = PAX*I_ESP_G4z_S_M3_vrr-PRX*I_ESP_G4z_S_M4_vrr;
        Double I_ESP_H5y_S_M3_vrr = PAY*I_ESP_G4y_S_M3_vrr-PRY*I_ESP_G4y_S_M4_vrr+4*oned2z*I_ESP_F3y_S_M3_vrr-4*oned2z*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_H4yz_S_M3_vrr = PAZ*I_ESP_G4y_S_M3_vrr-PRZ*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_H3y2z_S_M3_vrr = PAZ*I_ESP_G3yz_S_M3_vrr-PRZ*I_ESP_G3yz_S_M4_vrr+oned2z*I_ESP_F3y_S_M3_vrr-oned2z*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_Hy4z_S_M3_vrr = PAY*I_ESP_G4z_S_M3_vrr-PRY*I_ESP_G4z_S_M4_vrr;
        Double I_ESP_H5z_S_M3_vrr = PAZ*I_ESP_G4z_S_M3_vrr-PRZ*I_ESP_G4z_S_M4_vrr+4*oned2z*I_ESP_F3z_S_M3_vrr-4*oned2z*I_ESP_F3z_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 12 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S_M3
         * RHS shell quartet name: SQ_ESP_H_S_M4
         * RHS shell quartet name: SQ_ESP_G_S_M3
         * RHS shell quartet name: SQ_ESP_G_S_M4
         ************************************************************/
        Double I_ESP_I6x_S_M3_vrr = PAX*I_ESP_H5x_S_M3_vrr-PRX*I_ESP_H5x_S_M4_vrr+5*oned2z*I_ESP_G4x_S_M3_vrr-5*oned2z*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_I5xy_S_M3_vrr = PAY*I_ESP_H5x_S_M3_vrr-PRY*I_ESP_H5x_S_M4_vrr;
        Double I_ESP_I5xz_S_M3_vrr = PAZ*I_ESP_H5x_S_M3_vrr-PRZ*I_ESP_H5x_S_M4_vrr;
        Double I_ESP_I4x2y_S_M3_vrr = PAY*I_ESP_H4xy_S_M3_vrr-PRY*I_ESP_H4xy_S_M4_vrr+oned2z*I_ESP_G4x_S_M3_vrr-oned2z*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_I4x2z_S_M3_vrr = PAZ*I_ESP_H4xz_S_M3_vrr-PRZ*I_ESP_H4xz_S_M4_vrr+oned2z*I_ESP_G4x_S_M3_vrr-oned2z*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_I3x3y_S_M3_vrr = PAY*I_ESP_H3x2y_S_M3_vrr-PRY*I_ESP_H3x2y_S_M4_vrr+2*oned2z*I_ESP_G3xy_S_M3_vrr-2*oned2z*I_ESP_G3xy_S_M4_vrr;
        Double I_ESP_I2x4y_S_M3_vrr = PAX*I_ESP_Hx4y_S_M3_vrr-PRX*I_ESP_Hx4y_S_M4_vrr+oned2z*I_ESP_G4y_S_M3_vrr-oned2z*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_I2x4z_S_M3_vrr = PAX*I_ESP_Hx4z_S_M3_vrr-PRX*I_ESP_Hx4z_S_M4_vrr+oned2z*I_ESP_G4z_S_M3_vrr-oned2z*I_ESP_G4z_S_M4_vrr;
        Double I_ESP_Ix5y_S_M3_vrr = PAX*I_ESP_H5y_S_M3_vrr-PRX*I_ESP_H5y_S_M4_vrr;
        Double I_ESP_Ix5z_S_M3_vrr = PAX*I_ESP_H5z_S_M3_vrr-PRX*I_ESP_H5z_S_M4_vrr;
        Double I_ESP_I6y_S_M3_vrr = PAY*I_ESP_H5y_S_M3_vrr-PRY*I_ESP_H5y_S_M4_vrr+5*oned2z*I_ESP_G4y_S_M3_vrr-5*oned2z*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_I5yz_S_M3_vrr = PAZ*I_ESP_H5y_S_M3_vrr-PRZ*I_ESP_H5y_S_M4_vrr;
        Double I_ESP_I4y2z_S_M3_vrr = PAZ*I_ESP_H4yz_S_M3_vrr-PRZ*I_ESP_H4yz_S_M4_vrr+oned2z*I_ESP_G4y_S_M3_vrr-oned2z*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_I2y4z_S_M3_vrr = PAY*I_ESP_Hy4z_S_M3_vrr-PRY*I_ESP_Hy4z_S_M4_vrr+oned2z*I_ESP_G4z_S_M3_vrr-oned2z*I_ESP_G4z_S_M4_vrr;
        Double I_ESP_Iy5z_S_M3_vrr = PAY*I_ESP_H5z_S_M3_vrr-PRY*I_ESP_H5z_S_M4_vrr;
        Double I_ESP_I6z_S_M3_vrr = PAZ*I_ESP_H5z_S_M3_vrr-PRZ*I_ESP_H5z_S_M4_vrr+5*oned2z*I_ESP_G4z_S_M3_vrr-5*oned2z*I_ESP_G4z_S_M4_vrr;

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
         * shell quartet name: SQ_ESP_H_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M2
         * RHS shell quartet name: SQ_ESP_G_S_M3
         * RHS shell quartet name: SQ_ESP_F_S_M2
         * RHS shell quartet name: SQ_ESP_F_S_M3
         ************************************************************/
        Double I_ESP_H5x_S_M2_vrr = PAX*I_ESP_G4x_S_M2_vrr-PRX*I_ESP_G4x_S_M3_vrr+4*oned2z*I_ESP_F3x_S_M2_vrr-4*oned2z*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_H4xy_S_M2_vrr = PAY*I_ESP_G4x_S_M2_vrr-PRY*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_H4xz_S_M2_vrr = PAZ*I_ESP_G4x_S_M2_vrr-PRZ*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_H3x2y_S_M2_vrr = PAY*I_ESP_G3xy_S_M2_vrr-PRY*I_ESP_G3xy_S_M3_vrr+oned2z*I_ESP_F3x_S_M2_vrr-oned2z*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_H3x2z_S_M2_vrr = PAZ*I_ESP_G3xz_S_M2_vrr-PRZ*I_ESP_G3xz_S_M3_vrr+oned2z*I_ESP_F3x_S_M2_vrr-oned2z*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_H2x3y_S_M2_vrr = PAX*I_ESP_Gx3y_S_M2_vrr-PRX*I_ESP_Gx3y_S_M3_vrr+oned2z*I_ESP_F3y_S_M2_vrr-oned2z*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_H2x3z_S_M2_vrr = PAX*I_ESP_Gx3z_S_M2_vrr-PRX*I_ESP_Gx3z_S_M3_vrr+oned2z*I_ESP_F3z_S_M2_vrr-oned2z*I_ESP_F3z_S_M3_vrr;
        Double I_ESP_Hx4y_S_M2_vrr = PAX*I_ESP_G4y_S_M2_vrr-PRX*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_Hx4z_S_M2_vrr = PAX*I_ESP_G4z_S_M2_vrr-PRX*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_H5y_S_M2_vrr = PAY*I_ESP_G4y_S_M2_vrr-PRY*I_ESP_G4y_S_M3_vrr+4*oned2z*I_ESP_F3y_S_M2_vrr-4*oned2z*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_H4yz_S_M2_vrr = PAZ*I_ESP_G4y_S_M2_vrr-PRZ*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_H3y2z_S_M2_vrr = PAZ*I_ESP_G3yz_S_M2_vrr-PRZ*I_ESP_G3yz_S_M3_vrr+oned2z*I_ESP_F3y_S_M2_vrr-oned2z*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_Hy4z_S_M2_vrr = PAY*I_ESP_G4z_S_M2_vrr-PRY*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_H5z_S_M2_vrr = PAZ*I_ESP_G4z_S_M2_vrr-PRZ*I_ESP_G4z_S_M3_vrr+4*oned2z*I_ESP_F3z_S_M2_vrr-4*oned2z*I_ESP_F3z_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 10 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S_M2
         * RHS shell quartet name: SQ_ESP_H_S_M3
         * RHS shell quartet name: SQ_ESP_G_S_M2
         * RHS shell quartet name: SQ_ESP_G_S_M3
         ************************************************************/
        Double I_ESP_I6x_S_M2_vrr = PAX*I_ESP_H5x_S_M2_vrr-PRX*I_ESP_H5x_S_M3_vrr+5*oned2z*I_ESP_G4x_S_M2_vrr-5*oned2z*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_I5xy_S_M2_vrr = PAY*I_ESP_H5x_S_M2_vrr-PRY*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_I5xz_S_M2_vrr = PAZ*I_ESP_H5x_S_M2_vrr-PRZ*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_I4x2y_S_M2_vrr = PAY*I_ESP_H4xy_S_M2_vrr-PRY*I_ESP_H4xy_S_M3_vrr+oned2z*I_ESP_G4x_S_M2_vrr-oned2z*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_I4x2z_S_M2_vrr = PAZ*I_ESP_H4xz_S_M2_vrr-PRZ*I_ESP_H4xz_S_M3_vrr+oned2z*I_ESP_G4x_S_M2_vrr-oned2z*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_I3x3y_S_M2_vrr = PAY*I_ESP_H3x2y_S_M2_vrr-PRY*I_ESP_H3x2y_S_M3_vrr+2*oned2z*I_ESP_G3xy_S_M2_vrr-2*oned2z*I_ESP_G3xy_S_M3_vrr;
        Double I_ESP_I3x3z_S_M2_vrr = PAZ*I_ESP_H3x2z_S_M2_vrr-PRZ*I_ESP_H3x2z_S_M3_vrr+2*oned2z*I_ESP_G3xz_S_M2_vrr-2*oned2z*I_ESP_G3xz_S_M3_vrr;
        Double I_ESP_I2x4y_S_M2_vrr = PAX*I_ESP_Hx4y_S_M2_vrr-PRX*I_ESP_Hx4y_S_M3_vrr+oned2z*I_ESP_G4y_S_M2_vrr-oned2z*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_I2x4z_S_M2_vrr = PAX*I_ESP_Hx4z_S_M2_vrr-PRX*I_ESP_Hx4z_S_M3_vrr+oned2z*I_ESP_G4z_S_M2_vrr-oned2z*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_Ix5y_S_M2_vrr = PAX*I_ESP_H5y_S_M2_vrr-PRX*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_Ix5z_S_M2_vrr = PAX*I_ESP_H5z_S_M2_vrr-PRX*I_ESP_H5z_S_M3_vrr;
        Double I_ESP_I6y_S_M2_vrr = PAY*I_ESP_H5y_S_M2_vrr-PRY*I_ESP_H5y_S_M3_vrr+5*oned2z*I_ESP_G4y_S_M2_vrr-5*oned2z*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_I5yz_S_M2_vrr = PAZ*I_ESP_H5y_S_M2_vrr-PRZ*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_I4y2z_S_M2_vrr = PAZ*I_ESP_H4yz_S_M2_vrr-PRZ*I_ESP_H4yz_S_M3_vrr+oned2z*I_ESP_G4y_S_M2_vrr-oned2z*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_I3y3z_S_M2_vrr = PAZ*I_ESP_H3y2z_S_M2_vrr-PRZ*I_ESP_H3y2z_S_M3_vrr+2*oned2z*I_ESP_G3yz_S_M2_vrr-2*oned2z*I_ESP_G3yz_S_M3_vrr;
        Double I_ESP_I2y4z_S_M2_vrr = PAY*I_ESP_Hy4z_S_M2_vrr-PRY*I_ESP_Hy4z_S_M3_vrr+oned2z*I_ESP_G4z_S_M2_vrr-oned2z*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_Iy5z_S_M2_vrr = PAY*I_ESP_H5z_S_M2_vrr-PRY*I_ESP_H5z_S_M3_vrr;
        Double I_ESP_I6z_S_M2_vrr = PAZ*I_ESP_H5z_S_M2_vrr-PRZ*I_ESP_H5z_S_M3_vrr+5*oned2z*I_ESP_G4z_S_M2_vrr-5*oned2z*I_ESP_G4z_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 14 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_I_S_M2
         * RHS shell quartet name: SQ_ESP_I_S_M3
         * RHS shell quartet name: SQ_ESP_H_S_M2
         * RHS shell quartet name: SQ_ESP_H_S_M3
         ************************************************************/
        Double I_ESP_K7x_S_M2_vrr = PAX*I_ESP_I6x_S_M2_vrr-PRX*I_ESP_I6x_S_M3_vrr+6*oned2z*I_ESP_H5x_S_M2_vrr-6*oned2z*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_K6xy_S_M2_vrr = PAY*I_ESP_I6x_S_M2_vrr-PRY*I_ESP_I6x_S_M3_vrr;
        Double I_ESP_K6xz_S_M2_vrr = PAZ*I_ESP_I6x_S_M2_vrr-PRZ*I_ESP_I6x_S_M3_vrr;
        Double I_ESP_K5x2y_S_M2_vrr = PAY*I_ESP_I5xy_S_M2_vrr-PRY*I_ESP_I5xy_S_M3_vrr+oned2z*I_ESP_H5x_S_M2_vrr-oned2z*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_K5x2z_S_M2_vrr = PAZ*I_ESP_I5xz_S_M2_vrr-PRZ*I_ESP_I5xz_S_M3_vrr+oned2z*I_ESP_H5x_S_M2_vrr-oned2z*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_K4x3y_S_M2_vrr = PAY*I_ESP_I4x2y_S_M2_vrr-PRY*I_ESP_I4x2y_S_M3_vrr+2*oned2z*I_ESP_H4xy_S_M2_vrr-2*oned2z*I_ESP_H4xy_S_M3_vrr;
        Double I_ESP_K4x3z_S_M2_vrr = PAZ*I_ESP_I4x2z_S_M2_vrr-PRZ*I_ESP_I4x2z_S_M3_vrr+2*oned2z*I_ESP_H4xz_S_M2_vrr-2*oned2z*I_ESP_H4xz_S_M3_vrr;
        Double I_ESP_K3x4y_S_M2_vrr = PAX*I_ESP_I2x4y_S_M2_vrr-PRX*I_ESP_I2x4y_S_M3_vrr+2*oned2z*I_ESP_Hx4y_S_M2_vrr-2*oned2z*I_ESP_Hx4y_S_M3_vrr;
        Double I_ESP_K3x3yz_S_M2_vrr = PAZ*I_ESP_I3x3y_S_M2_vrr-PRZ*I_ESP_I3x3y_S_M3_vrr;
        Double I_ESP_K3x4z_S_M2_vrr = PAX*I_ESP_I2x4z_S_M2_vrr-PRX*I_ESP_I2x4z_S_M3_vrr+2*oned2z*I_ESP_Hx4z_S_M2_vrr-2*oned2z*I_ESP_Hx4z_S_M3_vrr;
        Double I_ESP_K2x5y_S_M2_vrr = PAX*I_ESP_Ix5y_S_M2_vrr-PRX*I_ESP_Ix5y_S_M3_vrr+oned2z*I_ESP_H5y_S_M2_vrr-oned2z*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_K2x5z_S_M2_vrr = PAX*I_ESP_Ix5z_S_M2_vrr-PRX*I_ESP_Ix5z_S_M3_vrr+oned2z*I_ESP_H5z_S_M2_vrr-oned2z*I_ESP_H5z_S_M3_vrr;
        Double I_ESP_Kx6y_S_M2_vrr = PAX*I_ESP_I6y_S_M2_vrr-PRX*I_ESP_I6y_S_M3_vrr;
        Double I_ESP_Kx6z_S_M2_vrr = PAX*I_ESP_I6z_S_M2_vrr-PRX*I_ESP_I6z_S_M3_vrr;
        Double I_ESP_K7y_S_M2_vrr = PAY*I_ESP_I6y_S_M2_vrr-PRY*I_ESP_I6y_S_M3_vrr+6*oned2z*I_ESP_H5y_S_M2_vrr-6*oned2z*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_K6yz_S_M2_vrr = PAZ*I_ESP_I6y_S_M2_vrr-PRZ*I_ESP_I6y_S_M3_vrr;
        Double I_ESP_K5y2z_S_M2_vrr = PAZ*I_ESP_I5yz_S_M2_vrr-PRZ*I_ESP_I5yz_S_M3_vrr+oned2z*I_ESP_H5y_S_M2_vrr-oned2z*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_K4y3z_S_M2_vrr = PAZ*I_ESP_I4y2z_S_M2_vrr-PRZ*I_ESP_I4y2z_S_M3_vrr+2*oned2z*I_ESP_H4yz_S_M2_vrr-2*oned2z*I_ESP_H4yz_S_M3_vrr;
        Double I_ESP_K3y4z_S_M2_vrr = PAY*I_ESP_I2y4z_S_M2_vrr-PRY*I_ESP_I2y4z_S_M3_vrr+2*oned2z*I_ESP_Hy4z_S_M2_vrr-2*oned2z*I_ESP_Hy4z_S_M3_vrr;
        Double I_ESP_K2y5z_S_M2_vrr = PAY*I_ESP_Iy5z_S_M2_vrr-PRY*I_ESP_Iy5z_S_M3_vrr+oned2z*I_ESP_H5z_S_M2_vrr-oned2z*I_ESP_H5z_S_M3_vrr;
        Double I_ESP_Ky6z_S_M2_vrr = PAY*I_ESP_I6z_S_M2_vrr-PRY*I_ESP_I6z_S_M3_vrr;
        Double I_ESP_K7z_S_M2_vrr = PAZ*I_ESP_I6z_S_M2_vrr-PRZ*I_ESP_I6z_S_M3_vrr+6*oned2z*I_ESP_H5z_S_M2_vrr-6*oned2z*I_ESP_H5z_S_M3_vrr;

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
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_D2x_S_M1_vrr = PAX*I_ESP_Px_S_M1_vrr-PRX*I_ESP_Px_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
        Double I_ESP_D2y_S_M1_vrr = PAY*I_ESP_Py_S_M1_vrr-PRY*I_ESP_Py_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
        Double I_ESP_D2z_S_M1_vrr = PAZ*I_ESP_Pz_S_M1_vrr-PRZ*I_ESP_Pz_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 4 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         ************************************************************/
        Double I_ESP_F3x_S_M1_vrr = PAX*I_ESP_D2x_S_M1_vrr-PRX*I_ESP_D2x_S_M2_vrr+2*oned2z*I_ESP_Px_S_M1_vrr-2*oned2z*I_ESP_Px_S_M2_vrr;
        Double I_ESP_F2xy_S_M1_vrr = PAY*I_ESP_D2x_S_M1_vrr-PRY*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_F2xz_S_M1_vrr = PAZ*I_ESP_D2x_S_M1_vrr-PRZ*I_ESP_D2x_S_M2_vrr;
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
         * shell quartet name: SQ_ESP_I_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S_M1
         * RHS shell quartet name: SQ_ESP_H_S_M2
         * RHS shell quartet name: SQ_ESP_G_S_M1
         * RHS shell quartet name: SQ_ESP_G_S_M2
         ************************************************************/
        Double I_ESP_I6x_S_M1_vrr = PAX*I_ESP_H5x_S_M1_vrr-PRX*I_ESP_H5x_S_M2_vrr+5*oned2z*I_ESP_G4x_S_M1_vrr-5*oned2z*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_I5xy_S_M1_vrr = PAY*I_ESP_H5x_S_M1_vrr-PRY*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_I5xz_S_M1_vrr = PAZ*I_ESP_H5x_S_M1_vrr-PRZ*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_I4x2y_S_M1_vrr = PAY*I_ESP_H4xy_S_M1_vrr-PRY*I_ESP_H4xy_S_M2_vrr+oned2z*I_ESP_G4x_S_M1_vrr-oned2z*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_I4x2z_S_M1_vrr = PAZ*I_ESP_H4xz_S_M1_vrr-PRZ*I_ESP_H4xz_S_M2_vrr+oned2z*I_ESP_G4x_S_M1_vrr-oned2z*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_I3x3y_S_M1_vrr = PAY*I_ESP_H3x2y_S_M1_vrr-PRY*I_ESP_H3x2y_S_M2_vrr+2*oned2z*I_ESP_G3xy_S_M1_vrr-2*oned2z*I_ESP_G3xy_S_M2_vrr;
        Double I_ESP_I3x2yz_S_M1_vrr = PAZ*I_ESP_H3x2y_S_M1_vrr-PRZ*I_ESP_H3x2y_S_M2_vrr;
        Double I_ESP_I3x3z_S_M1_vrr = PAZ*I_ESP_H3x2z_S_M1_vrr-PRZ*I_ESP_H3x2z_S_M2_vrr+2*oned2z*I_ESP_G3xz_S_M1_vrr-2*oned2z*I_ESP_G3xz_S_M2_vrr;
        Double I_ESP_I2x4y_S_M1_vrr = PAX*I_ESP_Hx4y_S_M1_vrr-PRX*I_ESP_Hx4y_S_M2_vrr+oned2z*I_ESP_G4y_S_M1_vrr-oned2z*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_I2x3yz_S_M1_vrr = PAZ*I_ESP_H2x3y_S_M1_vrr-PRZ*I_ESP_H2x3y_S_M2_vrr;
        Double I_ESP_I2xy3z_S_M1_vrr = PAY*I_ESP_H2x3z_S_M1_vrr-PRY*I_ESP_H2x3z_S_M2_vrr;
        Double I_ESP_I2x4z_S_M1_vrr = PAX*I_ESP_Hx4z_S_M1_vrr-PRX*I_ESP_Hx4z_S_M2_vrr+oned2z*I_ESP_G4z_S_M1_vrr-oned2z*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_Ix5y_S_M1_vrr = PAX*I_ESP_H5y_S_M1_vrr-PRX*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_Ix5z_S_M1_vrr = PAX*I_ESP_H5z_S_M1_vrr-PRX*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_I6y_S_M1_vrr = PAY*I_ESP_H5y_S_M1_vrr-PRY*I_ESP_H5y_S_M2_vrr+5*oned2z*I_ESP_G4y_S_M1_vrr-5*oned2z*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_I5yz_S_M1_vrr = PAZ*I_ESP_H5y_S_M1_vrr-PRZ*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_I4y2z_S_M1_vrr = PAZ*I_ESP_H4yz_S_M1_vrr-PRZ*I_ESP_H4yz_S_M2_vrr+oned2z*I_ESP_G4y_S_M1_vrr-oned2z*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_I3y3z_S_M1_vrr = PAZ*I_ESP_H3y2z_S_M1_vrr-PRZ*I_ESP_H3y2z_S_M2_vrr+2*oned2z*I_ESP_G3yz_S_M1_vrr-2*oned2z*I_ESP_G3yz_S_M2_vrr;
        Double I_ESP_I2y4z_S_M1_vrr = PAY*I_ESP_Hy4z_S_M1_vrr-PRY*I_ESP_Hy4z_S_M2_vrr+oned2z*I_ESP_G4z_S_M1_vrr-oned2z*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_Iy5z_S_M1_vrr = PAY*I_ESP_H5z_S_M1_vrr-PRY*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_I6z_S_M1_vrr = PAZ*I_ESP_H5z_S_M1_vrr-PRZ*I_ESP_H5z_S_M2_vrr+5*oned2z*I_ESP_G4z_S_M1_vrr-5*oned2z*I_ESP_G4z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 9 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_I_S_M1
         * RHS shell quartet name: SQ_ESP_I_S_M2
         * RHS shell quartet name: SQ_ESP_H_S_M1
         * RHS shell quartet name: SQ_ESP_H_S_M2
         ************************************************************/
        Double I_ESP_K7x_S_M1_vrr = PAX*I_ESP_I6x_S_M1_vrr-PRX*I_ESP_I6x_S_M2_vrr+6*oned2z*I_ESP_H5x_S_M1_vrr-6*oned2z*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_K6xy_S_M1_vrr = PAY*I_ESP_I6x_S_M1_vrr-PRY*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_K6xz_S_M1_vrr = PAZ*I_ESP_I6x_S_M1_vrr-PRZ*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_K5x2y_S_M1_vrr = PAY*I_ESP_I5xy_S_M1_vrr-PRY*I_ESP_I5xy_S_M2_vrr+oned2z*I_ESP_H5x_S_M1_vrr-oned2z*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_K5x2z_S_M1_vrr = PAZ*I_ESP_I5xz_S_M1_vrr-PRZ*I_ESP_I5xz_S_M2_vrr+oned2z*I_ESP_H5x_S_M1_vrr-oned2z*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_K4x3y_S_M1_vrr = PAY*I_ESP_I4x2y_S_M1_vrr-PRY*I_ESP_I4x2y_S_M2_vrr+2*oned2z*I_ESP_H4xy_S_M1_vrr-2*oned2z*I_ESP_H4xy_S_M2_vrr;
        Double I_ESP_K4x2yz_S_M1_vrr = PAZ*I_ESP_I4x2y_S_M1_vrr-PRZ*I_ESP_I4x2y_S_M2_vrr;
        Double I_ESP_K4x3z_S_M1_vrr = PAZ*I_ESP_I4x2z_S_M1_vrr-PRZ*I_ESP_I4x2z_S_M2_vrr+2*oned2z*I_ESP_H4xz_S_M1_vrr-2*oned2z*I_ESP_H4xz_S_M2_vrr;
        Double I_ESP_K3x4y_S_M1_vrr = PAX*I_ESP_I2x4y_S_M1_vrr-PRX*I_ESP_I2x4y_S_M2_vrr+2*oned2z*I_ESP_Hx4y_S_M1_vrr-2*oned2z*I_ESP_Hx4y_S_M2_vrr;
        Double I_ESP_K3x3yz_S_M1_vrr = PAZ*I_ESP_I3x3y_S_M1_vrr-PRZ*I_ESP_I3x3y_S_M2_vrr;
        Double I_ESP_K3xy3z_S_M1_vrr = PAY*I_ESP_I3x3z_S_M1_vrr-PRY*I_ESP_I3x3z_S_M2_vrr;
        Double I_ESP_K3x4z_S_M1_vrr = PAX*I_ESP_I2x4z_S_M1_vrr-PRX*I_ESP_I2x4z_S_M2_vrr+2*oned2z*I_ESP_Hx4z_S_M1_vrr-2*oned2z*I_ESP_Hx4z_S_M2_vrr;
        Double I_ESP_K2x5y_S_M1_vrr = PAX*I_ESP_Ix5y_S_M1_vrr-PRX*I_ESP_Ix5y_S_M2_vrr+oned2z*I_ESP_H5y_S_M1_vrr-oned2z*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_K2x4yz_S_M1_vrr = PAZ*I_ESP_I2x4y_S_M1_vrr-PRZ*I_ESP_I2x4y_S_M2_vrr;
        Double I_ESP_K2xy4z_S_M1_vrr = PAY*I_ESP_I2x4z_S_M1_vrr-PRY*I_ESP_I2x4z_S_M2_vrr;
        Double I_ESP_K2x5z_S_M1_vrr = PAX*I_ESP_Ix5z_S_M1_vrr-PRX*I_ESP_Ix5z_S_M2_vrr+oned2z*I_ESP_H5z_S_M1_vrr-oned2z*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_Kx6y_S_M1_vrr = PAX*I_ESP_I6y_S_M1_vrr-PRX*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_Kx3y3z_S_M1_vrr = PAX*I_ESP_I3y3z_S_M1_vrr-PRX*I_ESP_I3y3z_S_M2_vrr;
        Double I_ESP_Kx6z_S_M1_vrr = PAX*I_ESP_I6z_S_M1_vrr-PRX*I_ESP_I6z_S_M2_vrr;
        Double I_ESP_K7y_S_M1_vrr = PAY*I_ESP_I6y_S_M1_vrr-PRY*I_ESP_I6y_S_M2_vrr+6*oned2z*I_ESP_H5y_S_M1_vrr-6*oned2z*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_K6yz_S_M1_vrr = PAZ*I_ESP_I6y_S_M1_vrr-PRZ*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_K5y2z_S_M1_vrr = PAZ*I_ESP_I5yz_S_M1_vrr-PRZ*I_ESP_I5yz_S_M2_vrr+oned2z*I_ESP_H5y_S_M1_vrr-oned2z*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_K4y3z_S_M1_vrr = PAZ*I_ESP_I4y2z_S_M1_vrr-PRZ*I_ESP_I4y2z_S_M2_vrr+2*oned2z*I_ESP_H4yz_S_M1_vrr-2*oned2z*I_ESP_H4yz_S_M2_vrr;
        Double I_ESP_K3y4z_S_M1_vrr = PAY*I_ESP_I2y4z_S_M1_vrr-PRY*I_ESP_I2y4z_S_M2_vrr+2*oned2z*I_ESP_Hy4z_S_M1_vrr-2*oned2z*I_ESP_Hy4z_S_M2_vrr;
        Double I_ESP_K2y5z_S_M1_vrr = PAY*I_ESP_Iy5z_S_M1_vrr-PRY*I_ESP_Iy5z_S_M2_vrr+oned2z*I_ESP_H5z_S_M1_vrr-oned2z*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_Ky6z_S_M1_vrr = PAY*I_ESP_I6z_S_M1_vrr-PRY*I_ESP_I6z_S_M2_vrr;
        Double I_ESP_K7z_S_M1_vrr = PAZ*I_ESP_I6z_S_M1_vrr-PRZ*I_ESP_I6z_S_M2_vrr+6*oned2z*I_ESP_H5z_S_M1_vrr-6*oned2z*I_ESP_H5z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_L_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 11 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_K_S_M1
         * RHS shell quartet name: SQ_ESP_K_S_M2
         * RHS shell quartet name: SQ_ESP_I_S_M1
         * RHS shell quartet name: SQ_ESP_I_S_M2
         ************************************************************/
        Double I_ESP_L8x_S_M1_vrr = PAX*I_ESP_K7x_S_M1_vrr-PRX*I_ESP_K7x_S_M2_vrr+7*oned2z*I_ESP_I6x_S_M1_vrr-7*oned2z*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_L7xy_S_M1_vrr = PAY*I_ESP_K7x_S_M1_vrr-PRY*I_ESP_K7x_S_M2_vrr;
        Double I_ESP_L7xz_S_M1_vrr = PAZ*I_ESP_K7x_S_M1_vrr-PRZ*I_ESP_K7x_S_M2_vrr;
        Double I_ESP_L6x2y_S_M1_vrr = PAY*I_ESP_K6xy_S_M1_vrr-PRY*I_ESP_K6xy_S_M2_vrr+oned2z*I_ESP_I6x_S_M1_vrr-oned2z*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_L6x2z_S_M1_vrr = PAZ*I_ESP_K6xz_S_M1_vrr-PRZ*I_ESP_K6xz_S_M2_vrr+oned2z*I_ESP_I6x_S_M1_vrr-oned2z*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_L5x3y_S_M1_vrr = PAY*I_ESP_K5x2y_S_M1_vrr-PRY*I_ESP_K5x2y_S_M2_vrr+2*oned2z*I_ESP_I5xy_S_M1_vrr-2*oned2z*I_ESP_I5xy_S_M2_vrr;
        Double I_ESP_L5x2yz_S_M1_vrr = PAZ*I_ESP_K5x2y_S_M1_vrr-PRZ*I_ESP_K5x2y_S_M2_vrr;
        Double I_ESP_L5x3z_S_M1_vrr = PAZ*I_ESP_K5x2z_S_M1_vrr-PRZ*I_ESP_K5x2z_S_M2_vrr+2*oned2z*I_ESP_I5xz_S_M1_vrr-2*oned2z*I_ESP_I5xz_S_M2_vrr;
        Double I_ESP_L4x4y_S_M1_vrr = PAY*I_ESP_K4x3y_S_M1_vrr-PRY*I_ESP_K4x3y_S_M2_vrr+3*oned2z*I_ESP_I4x2y_S_M1_vrr-3*oned2z*I_ESP_I4x2y_S_M2_vrr;
        Double I_ESP_L4x3yz_S_M1_vrr = PAZ*I_ESP_K4x3y_S_M1_vrr-PRZ*I_ESP_K4x3y_S_M2_vrr;
        Double I_ESP_L4xy3z_S_M1_vrr = PAY*I_ESP_K4x3z_S_M1_vrr-PRY*I_ESP_K4x3z_S_M2_vrr;
        Double I_ESP_L4x4z_S_M1_vrr = PAZ*I_ESP_K4x3z_S_M1_vrr-PRZ*I_ESP_K4x3z_S_M2_vrr+3*oned2z*I_ESP_I4x2z_S_M1_vrr-3*oned2z*I_ESP_I4x2z_S_M2_vrr;
        Double I_ESP_L3x5y_S_M1_vrr = PAX*I_ESP_K2x5y_S_M1_vrr-PRX*I_ESP_K2x5y_S_M2_vrr+2*oned2z*I_ESP_Ix5y_S_M1_vrr-2*oned2z*I_ESP_Ix5y_S_M2_vrr;
        Double I_ESP_L3x4yz_S_M1_vrr = PAZ*I_ESP_K3x4y_S_M1_vrr-PRZ*I_ESP_K3x4y_S_M2_vrr;
        Double I_ESP_L3x3y2z_S_M1_vrr = PAZ*I_ESP_K3x3yz_S_M1_vrr-PRZ*I_ESP_K3x3yz_S_M2_vrr+oned2z*I_ESP_I3x3y_S_M1_vrr-oned2z*I_ESP_I3x3y_S_M2_vrr;
        Double I_ESP_L3xy4z_S_M1_vrr = PAY*I_ESP_K3x4z_S_M1_vrr-PRY*I_ESP_K3x4z_S_M2_vrr;
        Double I_ESP_L3x5z_S_M1_vrr = PAX*I_ESP_K2x5z_S_M1_vrr-PRX*I_ESP_K2x5z_S_M2_vrr+2*oned2z*I_ESP_Ix5z_S_M1_vrr-2*oned2z*I_ESP_Ix5z_S_M2_vrr;
        Double I_ESP_L2x6y_S_M1_vrr = PAX*I_ESP_Kx6y_S_M1_vrr-PRX*I_ESP_Kx6y_S_M2_vrr+oned2z*I_ESP_I6y_S_M1_vrr-oned2z*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_L2x5yz_S_M1_vrr = PAZ*I_ESP_K2x5y_S_M1_vrr-PRZ*I_ESP_K2x5y_S_M2_vrr;
        Double I_ESP_L2xy5z_S_M1_vrr = PAY*I_ESP_K2x5z_S_M1_vrr-PRY*I_ESP_K2x5z_S_M2_vrr;
        Double I_ESP_L2x6z_S_M1_vrr = PAX*I_ESP_Kx6z_S_M1_vrr-PRX*I_ESP_Kx6z_S_M2_vrr+oned2z*I_ESP_I6z_S_M1_vrr-oned2z*I_ESP_I6z_S_M2_vrr;
        Double I_ESP_Lx7y_S_M1_vrr = PAX*I_ESP_K7y_S_M1_vrr-PRX*I_ESP_K7y_S_M2_vrr;
        Double I_ESP_Lx4y3z_S_M1_vrr = PAX*I_ESP_K4y3z_S_M1_vrr-PRX*I_ESP_K4y3z_S_M2_vrr;
        Double I_ESP_Lx3y4z_S_M1_vrr = PAX*I_ESP_K3y4z_S_M1_vrr-PRX*I_ESP_K3y4z_S_M2_vrr;
        Double I_ESP_Lx7z_S_M1_vrr = PAX*I_ESP_K7z_S_M1_vrr-PRX*I_ESP_K7z_S_M2_vrr;
        Double I_ESP_L8y_S_M1_vrr = PAY*I_ESP_K7y_S_M1_vrr-PRY*I_ESP_K7y_S_M2_vrr+7*oned2z*I_ESP_I6y_S_M1_vrr-7*oned2z*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_L7yz_S_M1_vrr = PAZ*I_ESP_K7y_S_M1_vrr-PRZ*I_ESP_K7y_S_M2_vrr;
        Double I_ESP_L6y2z_S_M1_vrr = PAZ*I_ESP_K6yz_S_M1_vrr-PRZ*I_ESP_K6yz_S_M2_vrr+oned2z*I_ESP_I6y_S_M1_vrr-oned2z*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_L5y3z_S_M1_vrr = PAZ*I_ESP_K5y2z_S_M1_vrr-PRZ*I_ESP_K5y2z_S_M2_vrr+2*oned2z*I_ESP_I5yz_S_M1_vrr-2*oned2z*I_ESP_I5yz_S_M2_vrr;
        Double I_ESP_L4y4z_S_M1_vrr = PAZ*I_ESP_K4y3z_S_M1_vrr-PRZ*I_ESP_K4y3z_S_M2_vrr+3*oned2z*I_ESP_I4y2z_S_M1_vrr-3*oned2z*I_ESP_I4y2z_S_M2_vrr;
        Double I_ESP_L3y5z_S_M1_vrr = PAY*I_ESP_K2y5z_S_M1_vrr-PRY*I_ESP_K2y5z_S_M2_vrr+2*oned2z*I_ESP_Iy5z_S_M1_vrr-2*oned2z*I_ESP_Iy5z_S_M2_vrr;
        Double I_ESP_L2y6z_S_M1_vrr = PAY*I_ESP_Ky6z_S_M1_vrr-PRY*I_ESP_Ky6z_S_M2_vrr+oned2z*I_ESP_I6z_S_M1_vrr-oned2z*I_ESP_I6z_S_M2_vrr;
        Double I_ESP_Ly7z_S_M1_vrr = PAY*I_ESP_K7z_S_M1_vrr-PRY*I_ESP_K7z_S_M2_vrr;
        Double I_ESP_L8z_S_M1_vrr = PAZ*I_ESP_K7z_S_M1_vrr-PRZ*I_ESP_K7z_S_M2_vrr+7*oned2z*I_ESP_I6z_S_M1_vrr-7*oned2z*I_ESP_I6z_S_M2_vrr;

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
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_D2x_S_vrr = PAX*I_ESP_Px_S_vrr-PRX*I_ESP_Px_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_D2y_S_vrr = PAY*I_ESP_Py_S_vrr-PRY*I_ESP_Py_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_D2z_S_vrr = PAZ*I_ESP_Pz_S_vrr-PRZ*I_ESP_Pz_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 4 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         ************************************************************/
        Double I_ESP_F3x_S_vrr = PAX*I_ESP_D2x_S_vrr-PRX*I_ESP_D2x_S_M1_vrr+2*oned2z*I_ESP_Px_S_vrr-2*oned2z*I_ESP_Px_S_M1_vrr;
        Double I_ESP_F2xy_S_vrr = PAY*I_ESP_D2x_S_vrr-PRY*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_F2xz_S_vrr = PAZ*I_ESP_D2x_S_vrr-PRZ*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_F3y_S_vrr = PAY*I_ESP_D2y_S_vrr-PRY*I_ESP_D2y_S_M1_vrr+2*oned2z*I_ESP_Py_S_vrr-2*oned2z*I_ESP_Py_S_M1_vrr;
        Double I_ESP_F2yz_S_vrr = PAZ*I_ESP_D2y_S_vrr-PRZ*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_F3z_S_vrr = PAZ*I_ESP_D2z_S_vrr-PRZ*I_ESP_D2z_S_M1_vrr+2*oned2z*I_ESP_Pz_S_vrr-2*oned2z*I_ESP_Pz_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         ************************************************************/
        Double I_ESP_G4x_S_vrr = PAX*I_ESP_F3x_S_vrr-PRX*I_ESP_F3x_S_M1_vrr+3*oned2z*I_ESP_D2x_S_vrr-3*oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G3xy_S_vrr = PAY*I_ESP_F3x_S_vrr-PRY*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G3xz_S_vrr = PAZ*I_ESP_F3x_S_vrr-PRZ*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G2x2y_S_vrr = PAY*I_ESP_F2xy_S_vrr-PRY*I_ESP_F2xy_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G2x2z_S_vrr = PAZ*I_ESP_F2xz_S_vrr-PRZ*I_ESP_F2xz_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Gx3y_S_vrr = PAX*I_ESP_F3y_S_vrr-PRX*I_ESP_F3y_S_M1_vrr;
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
         * shell quartet name: SQ_ESP_K_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_I_S
         * RHS shell quartet name: SQ_ESP_I_S_M1
         * RHS shell quartet name: SQ_ESP_H_S
         * RHS shell quartet name: SQ_ESP_H_S_M1
         ************************************************************/
        Double I_ESP_K7x_S_vrr = PAX*I_ESP_I6x_S_vrr-PRX*I_ESP_I6x_S_M1_vrr+6*oned2z*I_ESP_H5x_S_vrr-6*oned2z*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_K6xy_S_vrr = PAY*I_ESP_I6x_S_vrr-PRY*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_K6xz_S_vrr = PAZ*I_ESP_I6x_S_vrr-PRZ*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_K5x2y_S_vrr = PAY*I_ESP_I5xy_S_vrr-PRY*I_ESP_I5xy_S_M1_vrr+oned2z*I_ESP_H5x_S_vrr-oned2z*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_K5xyz_S_vrr = PAZ*I_ESP_I5xy_S_vrr-PRZ*I_ESP_I5xy_S_M1_vrr;
        Double I_ESP_K5x2z_S_vrr = PAZ*I_ESP_I5xz_S_vrr-PRZ*I_ESP_I5xz_S_M1_vrr+oned2z*I_ESP_H5x_S_vrr-oned2z*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_K4x3y_S_vrr = PAY*I_ESP_I4x2y_S_vrr-PRY*I_ESP_I4x2y_S_M1_vrr+2*oned2z*I_ESP_H4xy_S_vrr-2*oned2z*I_ESP_H4xy_S_M1_vrr;
        Double I_ESP_K4x2yz_S_vrr = PAZ*I_ESP_I4x2y_S_vrr-PRZ*I_ESP_I4x2y_S_M1_vrr;
        Double I_ESP_K4xy2z_S_vrr = PAY*I_ESP_I4x2z_S_vrr-PRY*I_ESP_I4x2z_S_M1_vrr;
        Double I_ESP_K4x3z_S_vrr = PAZ*I_ESP_I4x2z_S_vrr-PRZ*I_ESP_I4x2z_S_M1_vrr+2*oned2z*I_ESP_H4xz_S_vrr-2*oned2z*I_ESP_H4xz_S_M1_vrr;
        Double I_ESP_K3x4y_S_vrr = PAX*I_ESP_I2x4y_S_vrr-PRX*I_ESP_I2x4y_S_M1_vrr+2*oned2z*I_ESP_Hx4y_S_vrr-2*oned2z*I_ESP_Hx4y_S_M1_vrr;
        Double I_ESP_K3x3yz_S_vrr = PAZ*I_ESP_I3x3y_S_vrr-PRZ*I_ESP_I3x3y_S_M1_vrr;
        Double I_ESP_K3x2y2z_S_vrr = PAZ*I_ESP_I3x2yz_S_vrr-PRZ*I_ESP_I3x2yz_S_M1_vrr+oned2z*I_ESP_H3x2y_S_vrr-oned2z*I_ESP_H3x2y_S_M1_vrr;
        Double I_ESP_K3xy3z_S_vrr = PAY*I_ESP_I3x3z_S_vrr-PRY*I_ESP_I3x3z_S_M1_vrr;
        Double I_ESP_K3x4z_S_vrr = PAX*I_ESP_I2x4z_S_vrr-PRX*I_ESP_I2x4z_S_M1_vrr+2*oned2z*I_ESP_Hx4z_S_vrr-2*oned2z*I_ESP_Hx4z_S_M1_vrr;
        Double I_ESP_K2x5y_S_vrr = PAX*I_ESP_Ix5y_S_vrr-PRX*I_ESP_Ix5y_S_M1_vrr+oned2z*I_ESP_H5y_S_vrr-oned2z*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_K2x4yz_S_vrr = PAZ*I_ESP_I2x4y_S_vrr-PRZ*I_ESP_I2x4y_S_M1_vrr;
        Double I_ESP_K2x3y2z_S_vrr = PAZ*I_ESP_I2x3yz_S_vrr-PRZ*I_ESP_I2x3yz_S_M1_vrr+oned2z*I_ESP_H2x3y_S_vrr-oned2z*I_ESP_H2x3y_S_M1_vrr;
        Double I_ESP_K2x2y3z_S_vrr = PAY*I_ESP_I2xy3z_S_vrr-PRY*I_ESP_I2xy3z_S_M1_vrr+oned2z*I_ESP_H2x3z_S_vrr-oned2z*I_ESP_H2x3z_S_M1_vrr;
        Double I_ESP_K2xy4z_S_vrr = PAY*I_ESP_I2x4z_S_vrr-PRY*I_ESP_I2x4z_S_M1_vrr;
        Double I_ESP_K2x5z_S_vrr = PAX*I_ESP_Ix5z_S_vrr-PRX*I_ESP_Ix5z_S_M1_vrr+oned2z*I_ESP_H5z_S_vrr-oned2z*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_Kx6y_S_vrr = PAX*I_ESP_I6y_S_vrr-PRX*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_Kx5yz_S_vrr = PAZ*I_ESP_Ix5y_S_vrr-PRZ*I_ESP_Ix5y_S_M1_vrr;
        Double I_ESP_Kx4y2z_S_vrr = PAX*I_ESP_I4y2z_S_vrr-PRX*I_ESP_I4y2z_S_M1_vrr;
        Double I_ESP_Kx3y3z_S_vrr = PAX*I_ESP_I3y3z_S_vrr-PRX*I_ESP_I3y3z_S_M1_vrr;
        Double I_ESP_Kx2y4z_S_vrr = PAX*I_ESP_I2y4z_S_vrr-PRX*I_ESP_I2y4z_S_M1_vrr;
        Double I_ESP_Kxy5z_S_vrr = PAY*I_ESP_Ix5z_S_vrr-PRY*I_ESP_Ix5z_S_M1_vrr;
        Double I_ESP_Kx6z_S_vrr = PAX*I_ESP_I6z_S_vrr-PRX*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_K7y_S_vrr = PAY*I_ESP_I6y_S_vrr-PRY*I_ESP_I6y_S_M1_vrr+6*oned2z*I_ESP_H5y_S_vrr-6*oned2z*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_K6yz_S_vrr = PAZ*I_ESP_I6y_S_vrr-PRZ*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_K5y2z_S_vrr = PAZ*I_ESP_I5yz_S_vrr-PRZ*I_ESP_I5yz_S_M1_vrr+oned2z*I_ESP_H5y_S_vrr-oned2z*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_K4y3z_S_vrr = PAZ*I_ESP_I4y2z_S_vrr-PRZ*I_ESP_I4y2z_S_M1_vrr+2*oned2z*I_ESP_H4yz_S_vrr-2*oned2z*I_ESP_H4yz_S_M1_vrr;
        Double I_ESP_K3y4z_S_vrr = PAY*I_ESP_I2y4z_S_vrr-PRY*I_ESP_I2y4z_S_M1_vrr+2*oned2z*I_ESP_Hy4z_S_vrr-2*oned2z*I_ESP_Hy4z_S_M1_vrr;
        Double I_ESP_K2y5z_S_vrr = PAY*I_ESP_Iy5z_S_vrr-PRY*I_ESP_Iy5z_S_M1_vrr+oned2z*I_ESP_H5z_S_vrr-oned2z*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_Ky6z_S_vrr = PAY*I_ESP_I6z_S_vrr-PRY*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_K7z_S_vrr = PAZ*I_ESP_I6z_S_vrr-PRZ*I_ESP_I6z_S_M1_vrr+6*oned2z*I_ESP_H5z_S_vrr-6*oned2z*I_ESP_H5z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_L_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_K_S
         * RHS shell quartet name: SQ_ESP_K_S_M1
         * RHS shell quartet name: SQ_ESP_I_S
         * RHS shell quartet name: SQ_ESP_I_S_M1
         ************************************************************/
        Double I_ESP_L8x_S_vrr = PAX*I_ESP_K7x_S_vrr-PRX*I_ESP_K7x_S_M1_vrr+7*oned2z*I_ESP_I6x_S_vrr-7*oned2z*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_L7xy_S_vrr = PAY*I_ESP_K7x_S_vrr-PRY*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_L7xz_S_vrr = PAZ*I_ESP_K7x_S_vrr-PRZ*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_L6x2y_S_vrr = PAY*I_ESP_K6xy_S_vrr-PRY*I_ESP_K6xy_S_M1_vrr+oned2z*I_ESP_I6x_S_vrr-oned2z*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_L6xyz_S_vrr = PAZ*I_ESP_K6xy_S_vrr-PRZ*I_ESP_K6xy_S_M1_vrr;
        Double I_ESP_L6x2z_S_vrr = PAZ*I_ESP_K6xz_S_vrr-PRZ*I_ESP_K6xz_S_M1_vrr+oned2z*I_ESP_I6x_S_vrr-oned2z*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_L5x3y_S_vrr = PAY*I_ESP_K5x2y_S_vrr-PRY*I_ESP_K5x2y_S_M1_vrr+2*oned2z*I_ESP_I5xy_S_vrr-2*oned2z*I_ESP_I5xy_S_M1_vrr;
        Double I_ESP_L5x2yz_S_vrr = PAZ*I_ESP_K5x2y_S_vrr-PRZ*I_ESP_K5x2y_S_M1_vrr;
        Double I_ESP_L5xy2z_S_vrr = PAY*I_ESP_K5x2z_S_vrr-PRY*I_ESP_K5x2z_S_M1_vrr;
        Double I_ESP_L5x3z_S_vrr = PAZ*I_ESP_K5x2z_S_vrr-PRZ*I_ESP_K5x2z_S_M1_vrr+2*oned2z*I_ESP_I5xz_S_vrr-2*oned2z*I_ESP_I5xz_S_M1_vrr;
        Double I_ESP_L4x4y_S_vrr = PAY*I_ESP_K4x3y_S_vrr-PRY*I_ESP_K4x3y_S_M1_vrr+3*oned2z*I_ESP_I4x2y_S_vrr-3*oned2z*I_ESP_I4x2y_S_M1_vrr;
        Double I_ESP_L4x3yz_S_vrr = PAZ*I_ESP_K4x3y_S_vrr-PRZ*I_ESP_K4x3y_S_M1_vrr;
        Double I_ESP_L4x2y2z_S_vrr = PAZ*I_ESP_K4x2yz_S_vrr-PRZ*I_ESP_K4x2yz_S_M1_vrr+oned2z*I_ESP_I4x2y_S_vrr-oned2z*I_ESP_I4x2y_S_M1_vrr;
        Double I_ESP_L4xy3z_S_vrr = PAY*I_ESP_K4x3z_S_vrr-PRY*I_ESP_K4x3z_S_M1_vrr;
        Double I_ESP_L4x4z_S_vrr = PAZ*I_ESP_K4x3z_S_vrr-PRZ*I_ESP_K4x3z_S_M1_vrr+3*oned2z*I_ESP_I4x2z_S_vrr-3*oned2z*I_ESP_I4x2z_S_M1_vrr;
        Double I_ESP_L3x5y_S_vrr = PAX*I_ESP_K2x5y_S_vrr-PRX*I_ESP_K2x5y_S_M1_vrr+2*oned2z*I_ESP_Ix5y_S_vrr-2*oned2z*I_ESP_Ix5y_S_M1_vrr;
        Double I_ESP_L3x4yz_S_vrr = PAZ*I_ESP_K3x4y_S_vrr-PRZ*I_ESP_K3x4y_S_M1_vrr;
        Double I_ESP_L3x3y2z_S_vrr = PAZ*I_ESP_K3x3yz_S_vrr-PRZ*I_ESP_K3x3yz_S_M1_vrr+oned2z*I_ESP_I3x3y_S_vrr-oned2z*I_ESP_I3x3y_S_M1_vrr;
        Double I_ESP_L3x2y3z_S_vrr = PAY*I_ESP_K3xy3z_S_vrr-PRY*I_ESP_K3xy3z_S_M1_vrr+oned2z*I_ESP_I3x3z_S_vrr-oned2z*I_ESP_I3x3z_S_M1_vrr;
        Double I_ESP_L3xy4z_S_vrr = PAY*I_ESP_K3x4z_S_vrr-PRY*I_ESP_K3x4z_S_M1_vrr;
        Double I_ESP_L3x5z_S_vrr = PAX*I_ESP_K2x5z_S_vrr-PRX*I_ESP_K2x5z_S_M1_vrr+2*oned2z*I_ESP_Ix5z_S_vrr-2*oned2z*I_ESP_Ix5z_S_M1_vrr;
        Double I_ESP_L2x6y_S_vrr = PAX*I_ESP_Kx6y_S_vrr-PRX*I_ESP_Kx6y_S_M1_vrr+oned2z*I_ESP_I6y_S_vrr-oned2z*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_L2x5yz_S_vrr = PAZ*I_ESP_K2x5y_S_vrr-PRZ*I_ESP_K2x5y_S_M1_vrr;
        Double I_ESP_L2x4y2z_S_vrr = PAZ*I_ESP_K2x4yz_S_vrr-PRZ*I_ESP_K2x4yz_S_M1_vrr+oned2z*I_ESP_I2x4y_S_vrr-oned2z*I_ESP_I2x4y_S_M1_vrr;
        Double I_ESP_L2x3y3z_S_vrr = PAX*I_ESP_Kx3y3z_S_vrr-PRX*I_ESP_Kx3y3z_S_M1_vrr+oned2z*I_ESP_I3y3z_S_vrr-oned2z*I_ESP_I3y3z_S_M1_vrr;
        Double I_ESP_L2x2y4z_S_vrr = PAY*I_ESP_K2xy4z_S_vrr-PRY*I_ESP_K2xy4z_S_M1_vrr+oned2z*I_ESP_I2x4z_S_vrr-oned2z*I_ESP_I2x4z_S_M1_vrr;
        Double I_ESP_L2xy5z_S_vrr = PAY*I_ESP_K2x5z_S_vrr-PRY*I_ESP_K2x5z_S_M1_vrr;
        Double I_ESP_L2x6z_S_vrr = PAX*I_ESP_Kx6z_S_vrr-PRX*I_ESP_Kx6z_S_M1_vrr+oned2z*I_ESP_I6z_S_vrr-oned2z*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_Lx7y_S_vrr = PAX*I_ESP_K7y_S_vrr-PRX*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_Lx6yz_S_vrr = PAZ*I_ESP_Kx6y_S_vrr-PRZ*I_ESP_Kx6y_S_M1_vrr;
        Double I_ESP_Lx5y2z_S_vrr = PAX*I_ESP_K5y2z_S_vrr-PRX*I_ESP_K5y2z_S_M1_vrr;
        Double I_ESP_Lx4y3z_S_vrr = PAX*I_ESP_K4y3z_S_vrr-PRX*I_ESP_K4y3z_S_M1_vrr;
        Double I_ESP_Lx3y4z_S_vrr = PAX*I_ESP_K3y4z_S_vrr-PRX*I_ESP_K3y4z_S_M1_vrr;
        Double I_ESP_Lx2y5z_S_vrr = PAX*I_ESP_K2y5z_S_vrr-PRX*I_ESP_K2y5z_S_M1_vrr;
        Double I_ESP_Lxy6z_S_vrr = PAY*I_ESP_Kx6z_S_vrr-PRY*I_ESP_Kx6z_S_M1_vrr;
        Double I_ESP_Lx7z_S_vrr = PAX*I_ESP_K7z_S_vrr-PRX*I_ESP_K7z_S_M1_vrr;
        Double I_ESP_L8y_S_vrr = PAY*I_ESP_K7y_S_vrr-PRY*I_ESP_K7y_S_M1_vrr+7*oned2z*I_ESP_I6y_S_vrr-7*oned2z*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_L7yz_S_vrr = PAZ*I_ESP_K7y_S_vrr-PRZ*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_L6y2z_S_vrr = PAZ*I_ESP_K6yz_S_vrr-PRZ*I_ESP_K6yz_S_M1_vrr+oned2z*I_ESP_I6y_S_vrr-oned2z*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_L5y3z_S_vrr = PAZ*I_ESP_K5y2z_S_vrr-PRZ*I_ESP_K5y2z_S_M1_vrr+2*oned2z*I_ESP_I5yz_S_vrr-2*oned2z*I_ESP_I5yz_S_M1_vrr;
        Double I_ESP_L4y4z_S_vrr = PAZ*I_ESP_K4y3z_S_vrr-PRZ*I_ESP_K4y3z_S_M1_vrr+3*oned2z*I_ESP_I4y2z_S_vrr-3*oned2z*I_ESP_I4y2z_S_M1_vrr;
        Double I_ESP_L3y5z_S_vrr = PAY*I_ESP_K2y5z_S_vrr-PRY*I_ESP_K2y5z_S_M1_vrr+2*oned2z*I_ESP_Iy5z_S_vrr-2*oned2z*I_ESP_Iy5z_S_M1_vrr;
        Double I_ESP_L2y6z_S_vrr = PAY*I_ESP_Ky6z_S_vrr-PRY*I_ESP_Ky6z_S_M1_vrr+oned2z*I_ESP_I6z_S_vrr-oned2z*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_Ly7z_S_vrr = PAY*I_ESP_K7z_S_vrr-PRY*I_ESP_K7z_S_M1_vrr;
        Double I_ESP_L8z_S_vrr = PAZ*I_ESP_K7z_S_vrr-PRZ*I_ESP_K7z_S_M1_vrr+7*oned2z*I_ESP_I6z_S_vrr-7*oned2z*I_ESP_I6z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_M_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_L_S
         * RHS shell quartet name: SQ_ESP_L_S_M1
         * RHS shell quartet name: SQ_ESP_K_S
         * RHS shell quartet name: SQ_ESP_K_S_M1
         ************************************************************/
        Double I_ESP_M9x_S_vrr = PAX*I_ESP_L8x_S_vrr-PRX*I_ESP_L8x_S_M1_vrr+8*oned2z*I_ESP_K7x_S_vrr-8*oned2z*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_M8xy_S_vrr = PAY*I_ESP_L8x_S_vrr-PRY*I_ESP_L8x_S_M1_vrr;
        Double I_ESP_M8xz_S_vrr = PAZ*I_ESP_L8x_S_vrr-PRZ*I_ESP_L8x_S_M1_vrr;
        Double I_ESP_M7x2y_S_vrr = PAY*I_ESP_L7xy_S_vrr-PRY*I_ESP_L7xy_S_M1_vrr+oned2z*I_ESP_K7x_S_vrr-oned2z*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_M7xyz_S_vrr = PAZ*I_ESP_L7xy_S_vrr-PRZ*I_ESP_L7xy_S_M1_vrr;
        Double I_ESP_M7x2z_S_vrr = PAZ*I_ESP_L7xz_S_vrr-PRZ*I_ESP_L7xz_S_M1_vrr+oned2z*I_ESP_K7x_S_vrr-oned2z*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_M6x3y_S_vrr = PAY*I_ESP_L6x2y_S_vrr-PRY*I_ESP_L6x2y_S_M1_vrr+2*oned2z*I_ESP_K6xy_S_vrr-2*oned2z*I_ESP_K6xy_S_M1_vrr;
        Double I_ESP_M6x2yz_S_vrr = PAZ*I_ESP_L6x2y_S_vrr-PRZ*I_ESP_L6x2y_S_M1_vrr;
        Double I_ESP_M6xy2z_S_vrr = PAY*I_ESP_L6x2z_S_vrr-PRY*I_ESP_L6x2z_S_M1_vrr;
        Double I_ESP_M6x3z_S_vrr = PAZ*I_ESP_L6x2z_S_vrr-PRZ*I_ESP_L6x2z_S_M1_vrr+2*oned2z*I_ESP_K6xz_S_vrr-2*oned2z*I_ESP_K6xz_S_M1_vrr;
        Double I_ESP_M5x4y_S_vrr = PAY*I_ESP_L5x3y_S_vrr-PRY*I_ESP_L5x3y_S_M1_vrr+3*oned2z*I_ESP_K5x2y_S_vrr-3*oned2z*I_ESP_K5x2y_S_M1_vrr;
        Double I_ESP_M5x3yz_S_vrr = PAZ*I_ESP_L5x3y_S_vrr-PRZ*I_ESP_L5x3y_S_M1_vrr;
        Double I_ESP_M5x2y2z_S_vrr = PAZ*I_ESP_L5x2yz_S_vrr-PRZ*I_ESP_L5x2yz_S_M1_vrr+oned2z*I_ESP_K5x2y_S_vrr-oned2z*I_ESP_K5x2y_S_M1_vrr;
        Double I_ESP_M5xy3z_S_vrr = PAY*I_ESP_L5x3z_S_vrr-PRY*I_ESP_L5x3z_S_M1_vrr;
        Double I_ESP_M5x4z_S_vrr = PAZ*I_ESP_L5x3z_S_vrr-PRZ*I_ESP_L5x3z_S_M1_vrr+3*oned2z*I_ESP_K5x2z_S_vrr-3*oned2z*I_ESP_K5x2z_S_M1_vrr;
        Double I_ESP_M4x5y_S_vrr = PAX*I_ESP_L3x5y_S_vrr-PRX*I_ESP_L3x5y_S_M1_vrr+3*oned2z*I_ESP_K2x5y_S_vrr-3*oned2z*I_ESP_K2x5y_S_M1_vrr;
        Double I_ESP_M4x4yz_S_vrr = PAZ*I_ESP_L4x4y_S_vrr-PRZ*I_ESP_L4x4y_S_M1_vrr;
        Double I_ESP_M4x3y2z_S_vrr = PAZ*I_ESP_L4x3yz_S_vrr-PRZ*I_ESP_L4x3yz_S_M1_vrr+oned2z*I_ESP_K4x3y_S_vrr-oned2z*I_ESP_K4x3y_S_M1_vrr;
        Double I_ESP_M4x2y3z_S_vrr = PAY*I_ESP_L4xy3z_S_vrr-PRY*I_ESP_L4xy3z_S_M1_vrr+oned2z*I_ESP_K4x3z_S_vrr-oned2z*I_ESP_K4x3z_S_M1_vrr;
        Double I_ESP_M4xy4z_S_vrr = PAY*I_ESP_L4x4z_S_vrr-PRY*I_ESP_L4x4z_S_M1_vrr;
        Double I_ESP_M4x5z_S_vrr = PAX*I_ESP_L3x5z_S_vrr-PRX*I_ESP_L3x5z_S_M1_vrr+3*oned2z*I_ESP_K2x5z_S_vrr-3*oned2z*I_ESP_K2x5z_S_M1_vrr;
        Double I_ESP_M3x6y_S_vrr = PAX*I_ESP_L2x6y_S_vrr-PRX*I_ESP_L2x6y_S_M1_vrr+2*oned2z*I_ESP_Kx6y_S_vrr-2*oned2z*I_ESP_Kx6y_S_M1_vrr;
        Double I_ESP_M3x5yz_S_vrr = PAZ*I_ESP_L3x5y_S_vrr-PRZ*I_ESP_L3x5y_S_M1_vrr;
        Double I_ESP_M3x4y2z_S_vrr = PAZ*I_ESP_L3x4yz_S_vrr-PRZ*I_ESP_L3x4yz_S_M1_vrr+oned2z*I_ESP_K3x4y_S_vrr-oned2z*I_ESP_K3x4y_S_M1_vrr;
        Double I_ESP_M3x3y3z_S_vrr = PAZ*I_ESP_L3x3y2z_S_vrr-PRZ*I_ESP_L3x3y2z_S_M1_vrr+2*oned2z*I_ESP_K3x3yz_S_vrr-2*oned2z*I_ESP_K3x3yz_S_M1_vrr;
        Double I_ESP_M3x2y4z_S_vrr = PAY*I_ESP_L3xy4z_S_vrr-PRY*I_ESP_L3xy4z_S_M1_vrr+oned2z*I_ESP_K3x4z_S_vrr-oned2z*I_ESP_K3x4z_S_M1_vrr;
        Double I_ESP_M3xy5z_S_vrr = PAY*I_ESP_L3x5z_S_vrr-PRY*I_ESP_L3x5z_S_M1_vrr;
        Double I_ESP_M3x6z_S_vrr = PAX*I_ESP_L2x6z_S_vrr-PRX*I_ESP_L2x6z_S_M1_vrr+2*oned2z*I_ESP_Kx6z_S_vrr-2*oned2z*I_ESP_Kx6z_S_M1_vrr;
        Double I_ESP_M2x7y_S_vrr = PAX*I_ESP_Lx7y_S_vrr-PRX*I_ESP_Lx7y_S_M1_vrr+oned2z*I_ESP_K7y_S_vrr-oned2z*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_M2x6yz_S_vrr = PAZ*I_ESP_L2x6y_S_vrr-PRZ*I_ESP_L2x6y_S_M1_vrr;
        Double I_ESP_M2x5y2z_S_vrr = PAZ*I_ESP_L2x5yz_S_vrr-PRZ*I_ESP_L2x5yz_S_M1_vrr+oned2z*I_ESP_K2x5y_S_vrr-oned2z*I_ESP_K2x5y_S_M1_vrr;
        Double I_ESP_M2x4y3z_S_vrr = PAX*I_ESP_Lx4y3z_S_vrr-PRX*I_ESP_Lx4y3z_S_M1_vrr+oned2z*I_ESP_K4y3z_S_vrr-oned2z*I_ESP_K4y3z_S_M1_vrr;
        Double I_ESP_M2x3y4z_S_vrr = PAX*I_ESP_Lx3y4z_S_vrr-PRX*I_ESP_Lx3y4z_S_M1_vrr+oned2z*I_ESP_K3y4z_S_vrr-oned2z*I_ESP_K3y4z_S_M1_vrr;
        Double I_ESP_M2x2y5z_S_vrr = PAY*I_ESP_L2xy5z_S_vrr-PRY*I_ESP_L2xy5z_S_M1_vrr+oned2z*I_ESP_K2x5z_S_vrr-oned2z*I_ESP_K2x5z_S_M1_vrr;
        Double I_ESP_M2xy6z_S_vrr = PAY*I_ESP_L2x6z_S_vrr-PRY*I_ESP_L2x6z_S_M1_vrr;
        Double I_ESP_M2x7z_S_vrr = PAX*I_ESP_Lx7z_S_vrr-PRX*I_ESP_Lx7z_S_M1_vrr+oned2z*I_ESP_K7z_S_vrr-oned2z*I_ESP_K7z_S_M1_vrr;
        Double I_ESP_Mx8y_S_vrr = PAX*I_ESP_L8y_S_vrr-PRX*I_ESP_L8y_S_M1_vrr;
        Double I_ESP_Mx7yz_S_vrr = PAZ*I_ESP_Lx7y_S_vrr-PRZ*I_ESP_Lx7y_S_M1_vrr;
        Double I_ESP_Mx6y2z_S_vrr = PAX*I_ESP_L6y2z_S_vrr-PRX*I_ESP_L6y2z_S_M1_vrr;
        Double I_ESP_Mx5y3z_S_vrr = PAX*I_ESP_L5y3z_S_vrr-PRX*I_ESP_L5y3z_S_M1_vrr;
        Double I_ESP_Mx4y4z_S_vrr = PAX*I_ESP_L4y4z_S_vrr-PRX*I_ESP_L4y4z_S_M1_vrr;
        Double I_ESP_Mx3y5z_S_vrr = PAX*I_ESP_L3y5z_S_vrr-PRX*I_ESP_L3y5z_S_M1_vrr;
        Double I_ESP_Mx2y6z_S_vrr = PAX*I_ESP_L2y6z_S_vrr-PRX*I_ESP_L2y6z_S_M1_vrr;
        Double I_ESP_Mxy7z_S_vrr = PAY*I_ESP_Lx7z_S_vrr-PRY*I_ESP_Lx7z_S_M1_vrr;
        Double I_ESP_Mx8z_S_vrr = PAX*I_ESP_L8z_S_vrr-PRX*I_ESP_L8z_S_M1_vrr;
        Double I_ESP_M9y_S_vrr = PAY*I_ESP_L8y_S_vrr-PRY*I_ESP_L8y_S_M1_vrr+8*oned2z*I_ESP_K7y_S_vrr-8*oned2z*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_M8yz_S_vrr = PAZ*I_ESP_L8y_S_vrr-PRZ*I_ESP_L8y_S_M1_vrr;
        Double I_ESP_M7y2z_S_vrr = PAZ*I_ESP_L7yz_S_vrr-PRZ*I_ESP_L7yz_S_M1_vrr+oned2z*I_ESP_K7y_S_vrr-oned2z*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_M6y3z_S_vrr = PAZ*I_ESP_L6y2z_S_vrr-PRZ*I_ESP_L6y2z_S_M1_vrr+2*oned2z*I_ESP_K6yz_S_vrr-2*oned2z*I_ESP_K6yz_S_M1_vrr;
        Double I_ESP_M5y4z_S_vrr = PAZ*I_ESP_L5y3z_S_vrr-PRZ*I_ESP_L5y3z_S_M1_vrr+3*oned2z*I_ESP_K5y2z_S_vrr-3*oned2z*I_ESP_K5y2z_S_M1_vrr;
        Double I_ESP_M4y5z_S_vrr = PAY*I_ESP_L3y5z_S_vrr-PRY*I_ESP_L3y5z_S_M1_vrr+3*oned2z*I_ESP_K2y5z_S_vrr-3*oned2z*I_ESP_K2y5z_S_M1_vrr;
        Double I_ESP_M3y6z_S_vrr = PAY*I_ESP_L2y6z_S_vrr-PRY*I_ESP_L2y6z_S_M1_vrr+2*oned2z*I_ESP_Ky6z_S_vrr-2*oned2z*I_ESP_Ky6z_S_M1_vrr;
        Double I_ESP_M2y7z_S_vrr = PAY*I_ESP_Ly7z_S_vrr-PRY*I_ESP_Ly7z_S_M1_vrr+oned2z*I_ESP_K7z_S_vrr-oned2z*I_ESP_K7z_S_M1_vrr;
        Double I_ESP_My8z_S_vrr = PAY*I_ESP_L8z_S_vrr-PRY*I_ESP_L8z_S_M1_vrr;
        Double I_ESP_M9z_S_vrr = PAZ*I_ESP_L8z_S_vrr-PRZ*I_ESP_L8z_S_M1_vrr+8*oned2z*I_ESP_K7z_S_vrr-8*oned2z*I_ESP_K7z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_M_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_M9x_S += I_ESP_M9x_S_vrr;
        I_ESP_M8xy_S += I_ESP_M8xy_S_vrr;
        I_ESP_M8xz_S += I_ESP_M8xz_S_vrr;
        I_ESP_M7x2y_S += I_ESP_M7x2y_S_vrr;
        I_ESP_M7xyz_S += I_ESP_M7xyz_S_vrr;
        I_ESP_M7x2z_S += I_ESP_M7x2z_S_vrr;
        I_ESP_M6x3y_S += I_ESP_M6x3y_S_vrr;
        I_ESP_M6x2yz_S += I_ESP_M6x2yz_S_vrr;
        I_ESP_M6xy2z_S += I_ESP_M6xy2z_S_vrr;
        I_ESP_M6x3z_S += I_ESP_M6x3z_S_vrr;
        I_ESP_M5x4y_S += I_ESP_M5x4y_S_vrr;
        I_ESP_M5x3yz_S += I_ESP_M5x3yz_S_vrr;
        I_ESP_M5x2y2z_S += I_ESP_M5x2y2z_S_vrr;
        I_ESP_M5xy3z_S += I_ESP_M5xy3z_S_vrr;
        I_ESP_M5x4z_S += I_ESP_M5x4z_S_vrr;
        I_ESP_M4x5y_S += I_ESP_M4x5y_S_vrr;
        I_ESP_M4x4yz_S += I_ESP_M4x4yz_S_vrr;
        I_ESP_M4x3y2z_S += I_ESP_M4x3y2z_S_vrr;
        I_ESP_M4x2y3z_S += I_ESP_M4x2y3z_S_vrr;
        I_ESP_M4xy4z_S += I_ESP_M4xy4z_S_vrr;
        I_ESP_M4x5z_S += I_ESP_M4x5z_S_vrr;
        I_ESP_M3x6y_S += I_ESP_M3x6y_S_vrr;
        I_ESP_M3x5yz_S += I_ESP_M3x5yz_S_vrr;
        I_ESP_M3x4y2z_S += I_ESP_M3x4y2z_S_vrr;
        I_ESP_M3x3y3z_S += I_ESP_M3x3y3z_S_vrr;
        I_ESP_M3x2y4z_S += I_ESP_M3x2y4z_S_vrr;
        I_ESP_M3xy5z_S += I_ESP_M3xy5z_S_vrr;
        I_ESP_M3x6z_S += I_ESP_M3x6z_S_vrr;
        I_ESP_M2x7y_S += I_ESP_M2x7y_S_vrr;
        I_ESP_M2x6yz_S += I_ESP_M2x6yz_S_vrr;
        I_ESP_M2x5y2z_S += I_ESP_M2x5y2z_S_vrr;
        I_ESP_M2x4y3z_S += I_ESP_M2x4y3z_S_vrr;
        I_ESP_M2x3y4z_S += I_ESP_M2x3y4z_S_vrr;
        I_ESP_M2x2y5z_S += I_ESP_M2x2y5z_S_vrr;
        I_ESP_M2xy6z_S += I_ESP_M2xy6z_S_vrr;
        I_ESP_M2x7z_S += I_ESP_M2x7z_S_vrr;
        I_ESP_Mx8y_S += I_ESP_Mx8y_S_vrr;
        I_ESP_Mx7yz_S += I_ESP_Mx7yz_S_vrr;
        I_ESP_Mx6y2z_S += I_ESP_Mx6y2z_S_vrr;
        I_ESP_Mx5y3z_S += I_ESP_Mx5y3z_S_vrr;
        I_ESP_Mx4y4z_S += I_ESP_Mx4y4z_S_vrr;
        I_ESP_Mx3y5z_S += I_ESP_Mx3y5z_S_vrr;
        I_ESP_Mx2y6z_S += I_ESP_Mx2y6z_S_vrr;
        I_ESP_Mxy7z_S += I_ESP_Mxy7z_S_vrr;
        I_ESP_Mx8z_S += I_ESP_Mx8z_S_vrr;
        I_ESP_M9y_S += I_ESP_M9y_S_vrr;
        I_ESP_M8yz_S += I_ESP_M8yz_S_vrr;
        I_ESP_M7y2z_S += I_ESP_M7y2z_S_vrr;
        I_ESP_M6y3z_S += I_ESP_M6y3z_S_vrr;
        I_ESP_M5y4z_S += I_ESP_M5y4z_S_vrr;
        I_ESP_M4y5z_S += I_ESP_M4y5z_S_vrr;
        I_ESP_M3y6z_S += I_ESP_M3y6z_S_vrr;
        I_ESP_M2y7z_S += I_ESP_M2y7z_S_vrr;
        I_ESP_My8z_S += I_ESP_My8z_S_vrr;
        I_ESP_M9z_S += I_ESP_M9z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_L_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_L8x_S += I_ESP_L8x_S_vrr;
        I_ESP_L7xy_S += I_ESP_L7xy_S_vrr;
        I_ESP_L7xz_S += I_ESP_L7xz_S_vrr;
        I_ESP_L6x2y_S += I_ESP_L6x2y_S_vrr;
        I_ESP_L6xyz_S += I_ESP_L6xyz_S_vrr;
        I_ESP_L6x2z_S += I_ESP_L6x2z_S_vrr;
        I_ESP_L5x3y_S += I_ESP_L5x3y_S_vrr;
        I_ESP_L5x2yz_S += I_ESP_L5x2yz_S_vrr;
        I_ESP_L5xy2z_S += I_ESP_L5xy2z_S_vrr;
        I_ESP_L5x3z_S += I_ESP_L5x3z_S_vrr;
        I_ESP_L4x4y_S += I_ESP_L4x4y_S_vrr;
        I_ESP_L4x3yz_S += I_ESP_L4x3yz_S_vrr;
        I_ESP_L4x2y2z_S += I_ESP_L4x2y2z_S_vrr;
        I_ESP_L4xy3z_S += I_ESP_L4xy3z_S_vrr;
        I_ESP_L4x4z_S += I_ESP_L4x4z_S_vrr;
        I_ESP_L3x5y_S += I_ESP_L3x5y_S_vrr;
        I_ESP_L3x4yz_S += I_ESP_L3x4yz_S_vrr;
        I_ESP_L3x3y2z_S += I_ESP_L3x3y2z_S_vrr;
        I_ESP_L3x2y3z_S += I_ESP_L3x2y3z_S_vrr;
        I_ESP_L3xy4z_S += I_ESP_L3xy4z_S_vrr;
        I_ESP_L3x5z_S += I_ESP_L3x5z_S_vrr;
        I_ESP_L2x6y_S += I_ESP_L2x6y_S_vrr;
        I_ESP_L2x5yz_S += I_ESP_L2x5yz_S_vrr;
        I_ESP_L2x4y2z_S += I_ESP_L2x4y2z_S_vrr;
        I_ESP_L2x3y3z_S += I_ESP_L2x3y3z_S_vrr;
        I_ESP_L2x2y4z_S += I_ESP_L2x2y4z_S_vrr;
        I_ESP_L2xy5z_S += I_ESP_L2xy5z_S_vrr;
        I_ESP_L2x6z_S += I_ESP_L2x6z_S_vrr;
        I_ESP_Lx7y_S += I_ESP_Lx7y_S_vrr;
        I_ESP_Lx6yz_S += I_ESP_Lx6yz_S_vrr;
        I_ESP_Lx5y2z_S += I_ESP_Lx5y2z_S_vrr;
        I_ESP_Lx4y3z_S += I_ESP_Lx4y3z_S_vrr;
        I_ESP_Lx3y4z_S += I_ESP_Lx3y4z_S_vrr;
        I_ESP_Lx2y5z_S += I_ESP_Lx2y5z_S_vrr;
        I_ESP_Lxy6z_S += I_ESP_Lxy6z_S_vrr;
        I_ESP_Lx7z_S += I_ESP_Lx7z_S_vrr;
        I_ESP_L8y_S += I_ESP_L8y_S_vrr;
        I_ESP_L7yz_S += I_ESP_L7yz_S_vrr;
        I_ESP_L6y2z_S += I_ESP_L6y2z_S_vrr;
        I_ESP_L5y3z_S += I_ESP_L5y3z_S_vrr;
        I_ESP_L4y4z_S += I_ESP_L4y4z_S_vrr;
        I_ESP_L3y5z_S += I_ESP_L3y5z_S_vrr;
        I_ESP_L2y6z_S += I_ESP_L2y6z_S_vrr;
        I_ESP_Ly7z_S += I_ESP_Ly7z_S_vrr;
        I_ESP_L8z_S += I_ESP_L8z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_K7x_S += I_ESP_K7x_S_vrr;
        I_ESP_K6xy_S += I_ESP_K6xy_S_vrr;
        I_ESP_K6xz_S += I_ESP_K6xz_S_vrr;
        I_ESP_K5x2y_S += I_ESP_K5x2y_S_vrr;
        I_ESP_K5xyz_S += I_ESP_K5xyz_S_vrr;
        I_ESP_K5x2z_S += I_ESP_K5x2z_S_vrr;
        I_ESP_K4x3y_S += I_ESP_K4x3y_S_vrr;
        I_ESP_K4x2yz_S += I_ESP_K4x2yz_S_vrr;
        I_ESP_K4xy2z_S += I_ESP_K4xy2z_S_vrr;
        I_ESP_K4x3z_S += I_ESP_K4x3z_S_vrr;
        I_ESP_K3x4y_S += I_ESP_K3x4y_S_vrr;
        I_ESP_K3x3yz_S += I_ESP_K3x3yz_S_vrr;
        I_ESP_K3x2y2z_S += I_ESP_K3x2y2z_S_vrr;
        I_ESP_K3xy3z_S += I_ESP_K3xy3z_S_vrr;
        I_ESP_K3x4z_S += I_ESP_K3x4z_S_vrr;
        I_ESP_K2x5y_S += I_ESP_K2x5y_S_vrr;
        I_ESP_K2x4yz_S += I_ESP_K2x4yz_S_vrr;
        I_ESP_K2x3y2z_S += I_ESP_K2x3y2z_S_vrr;
        I_ESP_K2x2y3z_S += I_ESP_K2x2y3z_S_vrr;
        I_ESP_K2xy4z_S += I_ESP_K2xy4z_S_vrr;
        I_ESP_K2x5z_S += I_ESP_K2x5z_S_vrr;
        I_ESP_Kx6y_S += I_ESP_Kx6y_S_vrr;
        I_ESP_Kx5yz_S += I_ESP_Kx5yz_S_vrr;
        I_ESP_Kx4y2z_S += I_ESP_Kx4y2z_S_vrr;
        I_ESP_Kx3y3z_S += I_ESP_Kx3y3z_S_vrr;
        I_ESP_Kx2y4z_S += I_ESP_Kx2y4z_S_vrr;
        I_ESP_Kxy5z_S += I_ESP_Kxy5z_S_vrr;
        I_ESP_Kx6z_S += I_ESP_Kx6z_S_vrr;
        I_ESP_K7y_S += I_ESP_K7y_S_vrr;
        I_ESP_K6yz_S += I_ESP_K6yz_S_vrr;
        I_ESP_K5y2z_S += I_ESP_K5y2z_S_vrr;
        I_ESP_K4y3z_S += I_ESP_K4y3z_S_vrr;
        I_ESP_K3y4z_S += I_ESP_K3y4z_S_vrr;
        I_ESP_K2y5z_S += I_ESP_K2y5z_S_vrr;
        I_ESP_Ky6z_S += I_ESP_Ky6z_S_vrr;
        I_ESP_K7z_S += I_ESP_K7z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_I6x_S += I_ESP_I6x_S_vrr;
        I_ESP_I5xy_S += I_ESP_I5xy_S_vrr;
        I_ESP_I5xz_S += I_ESP_I5xz_S_vrr;
        I_ESP_I4x2y_S += I_ESP_I4x2y_S_vrr;
        I_ESP_I4xyz_S += I_ESP_I4xyz_S_vrr;
        I_ESP_I4x2z_S += I_ESP_I4x2z_S_vrr;
        I_ESP_I3x3y_S += I_ESP_I3x3y_S_vrr;
        I_ESP_I3x2yz_S += I_ESP_I3x2yz_S_vrr;
        I_ESP_I3xy2z_S += I_ESP_I3xy2z_S_vrr;
        I_ESP_I3x3z_S += I_ESP_I3x3z_S_vrr;
        I_ESP_I2x4y_S += I_ESP_I2x4y_S_vrr;
        I_ESP_I2x3yz_S += I_ESP_I2x3yz_S_vrr;
        I_ESP_I2x2y2z_S += I_ESP_I2x2y2z_S_vrr;
        I_ESP_I2xy3z_S += I_ESP_I2xy3z_S_vrr;
        I_ESP_I2x4z_S += I_ESP_I2x4z_S_vrr;
        I_ESP_Ix5y_S += I_ESP_Ix5y_S_vrr;
        I_ESP_Ix4yz_S += I_ESP_Ix4yz_S_vrr;
        I_ESP_Ix3y2z_S += I_ESP_Ix3y2z_S_vrr;
        I_ESP_Ix2y3z_S += I_ESP_Ix2y3z_S_vrr;
        I_ESP_Ixy4z_S += I_ESP_Ixy4z_S_vrr;
        I_ESP_Ix5z_S += I_ESP_Ix5z_S_vrr;
        I_ESP_I6y_S += I_ESP_I6y_S_vrr;
        I_ESP_I5yz_S += I_ESP_I5yz_S_vrr;
        I_ESP_I4y2z_S += I_ESP_I4y2z_S_vrr;
        I_ESP_I3y3z_S += I_ESP_I3y3z_S_vrr;
        I_ESP_I2y4z_S += I_ESP_I2y4z_S_vrr;
        I_ESP_Iy5z_S += I_ESP_Iy5z_S_vrr;
        I_ESP_I6z_S += I_ESP_I6z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_H5x_S += I_ESP_H5x_S_vrr;
        I_ESP_H4xy_S += I_ESP_H4xy_S_vrr;
        I_ESP_H4xz_S += I_ESP_H4xz_S_vrr;
        I_ESP_H3x2y_S += I_ESP_H3x2y_S_vrr;
        I_ESP_H3xyz_S += I_ESP_H3xyz_S_vrr;
        I_ESP_H3x2z_S += I_ESP_H3x2z_S_vrr;
        I_ESP_H2x3y_S += I_ESP_H2x3y_S_vrr;
        I_ESP_H2x2yz_S += I_ESP_H2x2yz_S_vrr;
        I_ESP_H2xy2z_S += I_ESP_H2xy2z_S_vrr;
        I_ESP_H2x3z_S += I_ESP_H2x3z_S_vrr;
        I_ESP_Hx4y_S += I_ESP_Hx4y_S_vrr;
        I_ESP_Hx3yz_S += I_ESP_Hx3yz_S_vrr;
        I_ESP_Hx2y2z_S += I_ESP_Hx2y2z_S_vrr;
        I_ESP_Hxy3z_S += I_ESP_Hxy3z_S_vrr;
        I_ESP_Hx4z_S += I_ESP_Hx4z_S_vrr;
        I_ESP_H5y_S += I_ESP_H5y_S_vrr;
        I_ESP_H4yz_S += I_ESP_H4yz_S_vrr;
        I_ESP_H3y2z_S += I_ESP_H3y2z_S_vrr;
        I_ESP_H2y3z_S += I_ESP_H2y3z_S_vrr;
        I_ESP_Hy4z_S += I_ESP_Hy4z_S_vrr;
        I_ESP_H5z_S += I_ESP_H5z_S_vrr;
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
     * shell quartet name: SQ_ESP_H_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S
     * RHS shell quartet name: SQ_ESP_H_S
     ************************************************************/
    Double I_ESP_H5x_Px = I_ESP_I6x_S+ABX*I_ESP_H5x_S;
    Double I_ESP_H4xy_Px = I_ESP_I5xy_S+ABX*I_ESP_H4xy_S;
    Double I_ESP_H4xz_Px = I_ESP_I5xz_S+ABX*I_ESP_H4xz_S;
    Double I_ESP_H3x2y_Px = I_ESP_I4x2y_S+ABX*I_ESP_H3x2y_S;
    Double I_ESP_H3xyz_Px = I_ESP_I4xyz_S+ABX*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Px = I_ESP_I4x2z_S+ABX*I_ESP_H3x2z_S;
    Double I_ESP_H2x3y_Px = I_ESP_I3x3y_S+ABX*I_ESP_H2x3y_S;
    Double I_ESP_H2x2yz_Px = I_ESP_I3x2yz_S+ABX*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Px = I_ESP_I3xy2z_S+ABX*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Px = I_ESP_I3x3z_S+ABX*I_ESP_H2x3z_S;
    Double I_ESP_Hx4y_Px = I_ESP_I2x4y_S+ABX*I_ESP_Hx4y_S;
    Double I_ESP_Hx3yz_Px = I_ESP_I2x3yz_S+ABX*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Px = I_ESP_I2x2y2z_S+ABX*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Px = I_ESP_I2xy3z_S+ABX*I_ESP_Hxy3z_S;
    Double I_ESP_Hx4z_Px = I_ESP_I2x4z_S+ABX*I_ESP_Hx4z_S;
    Double I_ESP_H5y_Px = I_ESP_Ix5y_S+ABX*I_ESP_H5y_S;
    Double I_ESP_H4yz_Px = I_ESP_Ix4yz_S+ABX*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Px = I_ESP_Ix3y2z_S+ABX*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Px = I_ESP_Ix2y3z_S+ABX*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Px = I_ESP_Ixy4z_S+ABX*I_ESP_Hy4z_S;
    Double I_ESP_H5z_Px = I_ESP_Ix5z_S+ABX*I_ESP_H5z_S;
    Double I_ESP_H5x_Py = I_ESP_I5xy_S+ABY*I_ESP_H5x_S;
    Double I_ESP_H4xy_Py = I_ESP_I4x2y_S+ABY*I_ESP_H4xy_S;
    Double I_ESP_H4xz_Py = I_ESP_I4xyz_S+ABY*I_ESP_H4xz_S;
    Double I_ESP_H3x2y_Py = I_ESP_I3x3y_S+ABY*I_ESP_H3x2y_S;
    Double I_ESP_H3xyz_Py = I_ESP_I3x2yz_S+ABY*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Py = I_ESP_I3xy2z_S+ABY*I_ESP_H3x2z_S;
    Double I_ESP_H2x3y_Py = I_ESP_I2x4y_S+ABY*I_ESP_H2x3y_S;
    Double I_ESP_H2x2yz_Py = I_ESP_I2x3yz_S+ABY*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Py = I_ESP_I2x2y2z_S+ABY*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Py = I_ESP_I2xy3z_S+ABY*I_ESP_H2x3z_S;
    Double I_ESP_Hx4y_Py = I_ESP_Ix5y_S+ABY*I_ESP_Hx4y_S;
    Double I_ESP_Hx3yz_Py = I_ESP_Ix4yz_S+ABY*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Py = I_ESP_Ix3y2z_S+ABY*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Py = I_ESP_Ix2y3z_S+ABY*I_ESP_Hxy3z_S;
    Double I_ESP_Hx4z_Py = I_ESP_Ixy4z_S+ABY*I_ESP_Hx4z_S;
    Double I_ESP_H5y_Py = I_ESP_I6y_S+ABY*I_ESP_H5y_S;
    Double I_ESP_H4yz_Py = I_ESP_I5yz_S+ABY*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Py = I_ESP_I4y2z_S+ABY*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Py = I_ESP_I3y3z_S+ABY*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Py = I_ESP_I2y4z_S+ABY*I_ESP_Hy4z_S;
    Double I_ESP_H5z_Py = I_ESP_Iy5z_S+ABY*I_ESP_H5z_S;
    Double I_ESP_H5x_Pz = I_ESP_I5xz_S+ABZ*I_ESP_H5x_S;
    Double I_ESP_H4xy_Pz = I_ESP_I4xyz_S+ABZ*I_ESP_H4xy_S;
    Double I_ESP_H4xz_Pz = I_ESP_I4x2z_S+ABZ*I_ESP_H4xz_S;
    Double I_ESP_H3x2y_Pz = I_ESP_I3x2yz_S+ABZ*I_ESP_H3x2y_S;
    Double I_ESP_H3xyz_Pz = I_ESP_I3xy2z_S+ABZ*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Pz = I_ESP_I3x3z_S+ABZ*I_ESP_H3x2z_S;
    Double I_ESP_H2x3y_Pz = I_ESP_I2x3yz_S+ABZ*I_ESP_H2x3y_S;
    Double I_ESP_H2x2yz_Pz = I_ESP_I2x2y2z_S+ABZ*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Pz = I_ESP_I2xy3z_S+ABZ*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Pz = I_ESP_I2x4z_S+ABZ*I_ESP_H2x3z_S;
    Double I_ESP_Hx4y_Pz = I_ESP_Ix4yz_S+ABZ*I_ESP_Hx4y_S;
    Double I_ESP_Hx3yz_Pz = I_ESP_Ix3y2z_S+ABZ*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Pz = I_ESP_Ix2y3z_S+ABZ*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Pz = I_ESP_Ixy4z_S+ABZ*I_ESP_Hxy3z_S;
    Double I_ESP_Hx4z_Pz = I_ESP_Ix5z_S+ABZ*I_ESP_Hx4z_S;
    Double I_ESP_H5y_Pz = I_ESP_I5yz_S+ABZ*I_ESP_H5y_S;
    Double I_ESP_H4yz_Pz = I_ESP_I4y2z_S+ABZ*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Pz = I_ESP_I3y3z_S+ABZ*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Pz = I_ESP_I2y4z_S+ABZ*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Pz = I_ESP_Iy5z_S+ABZ*I_ESP_Hy4z_S;
    Double I_ESP_H5z_Pz = I_ESP_I6z_S+ABZ*I_ESP_H5z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_I_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_S
     * RHS shell quartet name: SQ_ESP_I_S
     ************************************************************/
    Double I_ESP_I6x_Px = I_ESP_K7x_S+ABX*I_ESP_I6x_S;
    Double I_ESP_I5xy_Px = I_ESP_K6xy_S+ABX*I_ESP_I5xy_S;
    Double I_ESP_I5xz_Px = I_ESP_K6xz_S+ABX*I_ESP_I5xz_S;
    Double I_ESP_I4x2y_Px = I_ESP_K5x2y_S+ABX*I_ESP_I4x2y_S;
    Double I_ESP_I4xyz_Px = I_ESP_K5xyz_S+ABX*I_ESP_I4xyz_S;
    Double I_ESP_I4x2z_Px = I_ESP_K5x2z_S+ABX*I_ESP_I4x2z_S;
    Double I_ESP_I3x3y_Px = I_ESP_K4x3y_S+ABX*I_ESP_I3x3y_S;
    Double I_ESP_I3x2yz_Px = I_ESP_K4x2yz_S+ABX*I_ESP_I3x2yz_S;
    Double I_ESP_I3xy2z_Px = I_ESP_K4xy2z_S+ABX*I_ESP_I3xy2z_S;
    Double I_ESP_I3x3z_Px = I_ESP_K4x3z_S+ABX*I_ESP_I3x3z_S;
    Double I_ESP_I2x4y_Px = I_ESP_K3x4y_S+ABX*I_ESP_I2x4y_S;
    Double I_ESP_I2x3yz_Px = I_ESP_K3x3yz_S+ABX*I_ESP_I2x3yz_S;
    Double I_ESP_I2x2y2z_Px = I_ESP_K3x2y2z_S+ABX*I_ESP_I2x2y2z_S;
    Double I_ESP_I2xy3z_Px = I_ESP_K3xy3z_S+ABX*I_ESP_I2xy3z_S;
    Double I_ESP_I2x4z_Px = I_ESP_K3x4z_S+ABX*I_ESP_I2x4z_S;
    Double I_ESP_Ix5y_Px = I_ESP_K2x5y_S+ABX*I_ESP_Ix5y_S;
    Double I_ESP_Ix4yz_Px = I_ESP_K2x4yz_S+ABX*I_ESP_Ix4yz_S;
    Double I_ESP_Ix3y2z_Px = I_ESP_K2x3y2z_S+ABX*I_ESP_Ix3y2z_S;
    Double I_ESP_Ix2y3z_Px = I_ESP_K2x2y3z_S+ABX*I_ESP_Ix2y3z_S;
    Double I_ESP_Ixy4z_Px = I_ESP_K2xy4z_S+ABX*I_ESP_Ixy4z_S;
    Double I_ESP_Ix5z_Px = I_ESP_K2x5z_S+ABX*I_ESP_Ix5z_S;
    Double I_ESP_I6y_Px = I_ESP_Kx6y_S+ABX*I_ESP_I6y_S;
    Double I_ESP_I5yz_Px = I_ESP_Kx5yz_S+ABX*I_ESP_I5yz_S;
    Double I_ESP_I4y2z_Px = I_ESP_Kx4y2z_S+ABX*I_ESP_I4y2z_S;
    Double I_ESP_I3y3z_Px = I_ESP_Kx3y3z_S+ABX*I_ESP_I3y3z_S;
    Double I_ESP_I2y4z_Px = I_ESP_Kx2y4z_S+ABX*I_ESP_I2y4z_S;
    Double I_ESP_Iy5z_Px = I_ESP_Kxy5z_S+ABX*I_ESP_Iy5z_S;
    Double I_ESP_I6z_Px = I_ESP_Kx6z_S+ABX*I_ESP_I6z_S;
    Double I_ESP_I6x_Py = I_ESP_K6xy_S+ABY*I_ESP_I6x_S;
    Double I_ESP_I5xy_Py = I_ESP_K5x2y_S+ABY*I_ESP_I5xy_S;
    Double I_ESP_I5xz_Py = I_ESP_K5xyz_S+ABY*I_ESP_I5xz_S;
    Double I_ESP_I4x2y_Py = I_ESP_K4x3y_S+ABY*I_ESP_I4x2y_S;
    Double I_ESP_I4xyz_Py = I_ESP_K4x2yz_S+ABY*I_ESP_I4xyz_S;
    Double I_ESP_I4x2z_Py = I_ESP_K4xy2z_S+ABY*I_ESP_I4x2z_S;
    Double I_ESP_I3x3y_Py = I_ESP_K3x4y_S+ABY*I_ESP_I3x3y_S;
    Double I_ESP_I3x2yz_Py = I_ESP_K3x3yz_S+ABY*I_ESP_I3x2yz_S;
    Double I_ESP_I3xy2z_Py = I_ESP_K3x2y2z_S+ABY*I_ESP_I3xy2z_S;
    Double I_ESP_I3x3z_Py = I_ESP_K3xy3z_S+ABY*I_ESP_I3x3z_S;
    Double I_ESP_I2x4y_Py = I_ESP_K2x5y_S+ABY*I_ESP_I2x4y_S;
    Double I_ESP_I2x3yz_Py = I_ESP_K2x4yz_S+ABY*I_ESP_I2x3yz_S;
    Double I_ESP_I2x2y2z_Py = I_ESP_K2x3y2z_S+ABY*I_ESP_I2x2y2z_S;
    Double I_ESP_I2xy3z_Py = I_ESP_K2x2y3z_S+ABY*I_ESP_I2xy3z_S;
    Double I_ESP_I2x4z_Py = I_ESP_K2xy4z_S+ABY*I_ESP_I2x4z_S;
    Double I_ESP_Ix5y_Py = I_ESP_Kx6y_S+ABY*I_ESP_Ix5y_S;
    Double I_ESP_Ix4yz_Py = I_ESP_Kx5yz_S+ABY*I_ESP_Ix4yz_S;
    Double I_ESP_Ix3y2z_Py = I_ESP_Kx4y2z_S+ABY*I_ESP_Ix3y2z_S;
    Double I_ESP_Ix2y3z_Py = I_ESP_Kx3y3z_S+ABY*I_ESP_Ix2y3z_S;
    Double I_ESP_Ixy4z_Py = I_ESP_Kx2y4z_S+ABY*I_ESP_Ixy4z_S;
    Double I_ESP_Ix5z_Py = I_ESP_Kxy5z_S+ABY*I_ESP_Ix5z_S;
    Double I_ESP_I6y_Py = I_ESP_K7y_S+ABY*I_ESP_I6y_S;
    Double I_ESP_I5yz_Py = I_ESP_K6yz_S+ABY*I_ESP_I5yz_S;
    Double I_ESP_I4y2z_Py = I_ESP_K5y2z_S+ABY*I_ESP_I4y2z_S;
    Double I_ESP_I3y3z_Py = I_ESP_K4y3z_S+ABY*I_ESP_I3y3z_S;
    Double I_ESP_I2y4z_Py = I_ESP_K3y4z_S+ABY*I_ESP_I2y4z_S;
    Double I_ESP_Iy5z_Py = I_ESP_K2y5z_S+ABY*I_ESP_Iy5z_S;
    Double I_ESP_I6z_Py = I_ESP_Ky6z_S+ABY*I_ESP_I6z_S;
    Double I_ESP_I6x_Pz = I_ESP_K6xz_S+ABZ*I_ESP_I6x_S;
    Double I_ESP_I5xy_Pz = I_ESP_K5xyz_S+ABZ*I_ESP_I5xy_S;
    Double I_ESP_I5xz_Pz = I_ESP_K5x2z_S+ABZ*I_ESP_I5xz_S;
    Double I_ESP_I4x2y_Pz = I_ESP_K4x2yz_S+ABZ*I_ESP_I4x2y_S;
    Double I_ESP_I4xyz_Pz = I_ESP_K4xy2z_S+ABZ*I_ESP_I4xyz_S;
    Double I_ESP_I4x2z_Pz = I_ESP_K4x3z_S+ABZ*I_ESP_I4x2z_S;
    Double I_ESP_I3x3y_Pz = I_ESP_K3x3yz_S+ABZ*I_ESP_I3x3y_S;
    Double I_ESP_I3x2yz_Pz = I_ESP_K3x2y2z_S+ABZ*I_ESP_I3x2yz_S;
    Double I_ESP_I3xy2z_Pz = I_ESP_K3xy3z_S+ABZ*I_ESP_I3xy2z_S;
    Double I_ESP_I3x3z_Pz = I_ESP_K3x4z_S+ABZ*I_ESP_I3x3z_S;
    Double I_ESP_I2x4y_Pz = I_ESP_K2x4yz_S+ABZ*I_ESP_I2x4y_S;
    Double I_ESP_I2x3yz_Pz = I_ESP_K2x3y2z_S+ABZ*I_ESP_I2x3yz_S;
    Double I_ESP_I2x2y2z_Pz = I_ESP_K2x2y3z_S+ABZ*I_ESP_I2x2y2z_S;
    Double I_ESP_I2xy3z_Pz = I_ESP_K2xy4z_S+ABZ*I_ESP_I2xy3z_S;
    Double I_ESP_I2x4z_Pz = I_ESP_K2x5z_S+ABZ*I_ESP_I2x4z_S;
    Double I_ESP_Ix5y_Pz = I_ESP_Kx5yz_S+ABZ*I_ESP_Ix5y_S;
    Double I_ESP_Ix4yz_Pz = I_ESP_Kx4y2z_S+ABZ*I_ESP_Ix4yz_S;
    Double I_ESP_Ix3y2z_Pz = I_ESP_Kx3y3z_S+ABZ*I_ESP_Ix3y2z_S;
    Double I_ESP_Ix2y3z_Pz = I_ESP_Kx2y4z_S+ABZ*I_ESP_Ix2y3z_S;
    Double I_ESP_Ixy4z_Pz = I_ESP_Kxy5z_S+ABZ*I_ESP_Ixy4z_S;
    Double I_ESP_Ix5z_Pz = I_ESP_Kx6z_S+ABZ*I_ESP_Ix5z_S;
    Double I_ESP_I6y_Pz = I_ESP_K6yz_S+ABZ*I_ESP_I6y_S;
    Double I_ESP_I5yz_Pz = I_ESP_K5y2z_S+ABZ*I_ESP_I5yz_S;
    Double I_ESP_I4y2z_Pz = I_ESP_K4y3z_S+ABZ*I_ESP_I4y2z_S;
    Double I_ESP_I3y3z_Pz = I_ESP_K3y4z_S+ABZ*I_ESP_I3y3z_S;
    Double I_ESP_I2y4z_Pz = I_ESP_K2y5z_S+ABZ*I_ESP_I2y4z_S;
    Double I_ESP_Iy5z_Pz = I_ESP_Ky6z_S+ABZ*I_ESP_Iy5z_S;
    Double I_ESP_I6z_Pz = I_ESP_K7z_S+ABZ*I_ESP_I6z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 63 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P
     * RHS shell quartet name: SQ_ESP_H_P
     ************************************************************/
    Double I_ESP_H5x_D2x = I_ESP_I6x_Px+ABX*I_ESP_H5x_Px;
    Double I_ESP_H4xy_D2x = I_ESP_I5xy_Px+ABX*I_ESP_H4xy_Px;
    Double I_ESP_H4xz_D2x = I_ESP_I5xz_Px+ABX*I_ESP_H4xz_Px;
    Double I_ESP_H3x2y_D2x = I_ESP_I4x2y_Px+ABX*I_ESP_H3x2y_Px;
    Double I_ESP_H3xyz_D2x = I_ESP_I4xyz_Px+ABX*I_ESP_H3xyz_Px;
    Double I_ESP_H3x2z_D2x = I_ESP_I4x2z_Px+ABX*I_ESP_H3x2z_Px;
    Double I_ESP_H2x3y_D2x = I_ESP_I3x3y_Px+ABX*I_ESP_H2x3y_Px;
    Double I_ESP_H2x2yz_D2x = I_ESP_I3x2yz_Px+ABX*I_ESP_H2x2yz_Px;
    Double I_ESP_H2xy2z_D2x = I_ESP_I3xy2z_Px+ABX*I_ESP_H2xy2z_Px;
    Double I_ESP_H2x3z_D2x = I_ESP_I3x3z_Px+ABX*I_ESP_H2x3z_Px;
    Double I_ESP_Hx4y_D2x = I_ESP_I2x4y_Px+ABX*I_ESP_Hx4y_Px;
    Double I_ESP_Hx3yz_D2x = I_ESP_I2x3yz_Px+ABX*I_ESP_Hx3yz_Px;
    Double I_ESP_Hx2y2z_D2x = I_ESP_I2x2y2z_Px+ABX*I_ESP_Hx2y2z_Px;
    Double I_ESP_Hxy3z_D2x = I_ESP_I2xy3z_Px+ABX*I_ESP_Hxy3z_Px;
    Double I_ESP_Hx4z_D2x = I_ESP_I2x4z_Px+ABX*I_ESP_Hx4z_Px;
    Double I_ESP_H5y_D2x = I_ESP_Ix5y_Px+ABX*I_ESP_H5y_Px;
    Double I_ESP_H4yz_D2x = I_ESP_Ix4yz_Px+ABX*I_ESP_H4yz_Px;
    Double I_ESP_H3y2z_D2x = I_ESP_Ix3y2z_Px+ABX*I_ESP_H3y2z_Px;
    Double I_ESP_H2y3z_D2x = I_ESP_Ix2y3z_Px+ABX*I_ESP_H2y3z_Px;
    Double I_ESP_Hy4z_D2x = I_ESP_Ixy4z_Px+ABX*I_ESP_Hy4z_Px;
    Double I_ESP_H5z_D2x = I_ESP_Ix5z_Px+ABX*I_ESP_H5z_Px;
    Double I_ESP_H5x_D2y = I_ESP_I5xy_Py+ABY*I_ESP_H5x_Py;
    Double I_ESP_H4xy_D2y = I_ESP_I4x2y_Py+ABY*I_ESP_H4xy_Py;
    Double I_ESP_H4xz_D2y = I_ESP_I4xyz_Py+ABY*I_ESP_H4xz_Py;
    Double I_ESP_H3x2y_D2y = I_ESP_I3x3y_Py+ABY*I_ESP_H3x2y_Py;
    Double I_ESP_H3xyz_D2y = I_ESP_I3x2yz_Py+ABY*I_ESP_H3xyz_Py;
    Double I_ESP_H3x2z_D2y = I_ESP_I3xy2z_Py+ABY*I_ESP_H3x2z_Py;
    Double I_ESP_H2x3y_D2y = I_ESP_I2x4y_Py+ABY*I_ESP_H2x3y_Py;
    Double I_ESP_H2x2yz_D2y = I_ESP_I2x3yz_Py+ABY*I_ESP_H2x2yz_Py;
    Double I_ESP_H2xy2z_D2y = I_ESP_I2x2y2z_Py+ABY*I_ESP_H2xy2z_Py;
    Double I_ESP_H2x3z_D2y = I_ESP_I2xy3z_Py+ABY*I_ESP_H2x3z_Py;
    Double I_ESP_Hx4y_D2y = I_ESP_Ix5y_Py+ABY*I_ESP_Hx4y_Py;
    Double I_ESP_Hx3yz_D2y = I_ESP_Ix4yz_Py+ABY*I_ESP_Hx3yz_Py;
    Double I_ESP_Hx2y2z_D2y = I_ESP_Ix3y2z_Py+ABY*I_ESP_Hx2y2z_Py;
    Double I_ESP_Hxy3z_D2y = I_ESP_Ix2y3z_Py+ABY*I_ESP_Hxy3z_Py;
    Double I_ESP_Hx4z_D2y = I_ESP_Ixy4z_Py+ABY*I_ESP_Hx4z_Py;
    Double I_ESP_H5y_D2y = I_ESP_I6y_Py+ABY*I_ESP_H5y_Py;
    Double I_ESP_H4yz_D2y = I_ESP_I5yz_Py+ABY*I_ESP_H4yz_Py;
    Double I_ESP_H3y2z_D2y = I_ESP_I4y2z_Py+ABY*I_ESP_H3y2z_Py;
    Double I_ESP_H2y3z_D2y = I_ESP_I3y3z_Py+ABY*I_ESP_H2y3z_Py;
    Double I_ESP_Hy4z_D2y = I_ESP_I2y4z_Py+ABY*I_ESP_Hy4z_Py;
    Double I_ESP_H5z_D2y = I_ESP_Iy5z_Py+ABY*I_ESP_H5z_Py;
    Double I_ESP_H5x_D2z = I_ESP_I5xz_Pz+ABZ*I_ESP_H5x_Pz;
    Double I_ESP_H4xy_D2z = I_ESP_I4xyz_Pz+ABZ*I_ESP_H4xy_Pz;
    Double I_ESP_H4xz_D2z = I_ESP_I4x2z_Pz+ABZ*I_ESP_H4xz_Pz;
    Double I_ESP_H3x2y_D2z = I_ESP_I3x2yz_Pz+ABZ*I_ESP_H3x2y_Pz;
    Double I_ESP_H3xyz_D2z = I_ESP_I3xy2z_Pz+ABZ*I_ESP_H3xyz_Pz;
    Double I_ESP_H3x2z_D2z = I_ESP_I3x3z_Pz+ABZ*I_ESP_H3x2z_Pz;
    Double I_ESP_H2x3y_D2z = I_ESP_I2x3yz_Pz+ABZ*I_ESP_H2x3y_Pz;
    Double I_ESP_H2x2yz_D2z = I_ESP_I2x2y2z_Pz+ABZ*I_ESP_H2x2yz_Pz;
    Double I_ESP_H2xy2z_D2z = I_ESP_I2xy3z_Pz+ABZ*I_ESP_H2xy2z_Pz;
    Double I_ESP_H2x3z_D2z = I_ESP_I2x4z_Pz+ABZ*I_ESP_H2x3z_Pz;
    Double I_ESP_Hx4y_D2z = I_ESP_Ix4yz_Pz+ABZ*I_ESP_Hx4y_Pz;
    Double I_ESP_Hx3yz_D2z = I_ESP_Ix3y2z_Pz+ABZ*I_ESP_Hx3yz_Pz;
    Double I_ESP_Hx2y2z_D2z = I_ESP_Ix2y3z_Pz+ABZ*I_ESP_Hx2y2z_Pz;
    Double I_ESP_Hxy3z_D2z = I_ESP_Ixy4z_Pz+ABZ*I_ESP_Hxy3z_Pz;
    Double I_ESP_Hx4z_D2z = I_ESP_Ix5z_Pz+ABZ*I_ESP_Hx4z_Pz;
    Double I_ESP_H5y_D2z = I_ESP_I5yz_Pz+ABZ*I_ESP_H5y_Pz;
    Double I_ESP_H4yz_D2z = I_ESP_I4y2z_Pz+ABZ*I_ESP_H4yz_Pz;
    Double I_ESP_H3y2z_D2z = I_ESP_I3y3z_Pz+ABZ*I_ESP_H3y2z_Pz;
    Double I_ESP_H2y3z_D2z = I_ESP_I2y4z_Pz+ABZ*I_ESP_H2y3z_Pz;
    Double I_ESP_Hy4z_D2z = I_ESP_Iy5z_Pz+ABZ*I_ESP_Hy4z_Pz;
    Double I_ESP_H5z_D2z = I_ESP_I6z_Pz+ABZ*I_ESP_H5z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_K_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_S
     * RHS shell quartet name: SQ_ESP_K_S
     ************************************************************/
    Double I_ESP_K7x_Px = I_ESP_L8x_S+ABX*I_ESP_K7x_S;
    Double I_ESP_K6xy_Px = I_ESP_L7xy_S+ABX*I_ESP_K6xy_S;
    Double I_ESP_K6xz_Px = I_ESP_L7xz_S+ABX*I_ESP_K6xz_S;
    Double I_ESP_K5x2y_Px = I_ESP_L6x2y_S+ABX*I_ESP_K5x2y_S;
    Double I_ESP_K5xyz_Px = I_ESP_L6xyz_S+ABX*I_ESP_K5xyz_S;
    Double I_ESP_K5x2z_Px = I_ESP_L6x2z_S+ABX*I_ESP_K5x2z_S;
    Double I_ESP_K4x3y_Px = I_ESP_L5x3y_S+ABX*I_ESP_K4x3y_S;
    Double I_ESP_K4x2yz_Px = I_ESP_L5x2yz_S+ABX*I_ESP_K4x2yz_S;
    Double I_ESP_K4xy2z_Px = I_ESP_L5xy2z_S+ABX*I_ESP_K4xy2z_S;
    Double I_ESP_K4x3z_Px = I_ESP_L5x3z_S+ABX*I_ESP_K4x3z_S;
    Double I_ESP_K3x4y_Px = I_ESP_L4x4y_S+ABX*I_ESP_K3x4y_S;
    Double I_ESP_K3x3yz_Px = I_ESP_L4x3yz_S+ABX*I_ESP_K3x3yz_S;
    Double I_ESP_K3x2y2z_Px = I_ESP_L4x2y2z_S+ABX*I_ESP_K3x2y2z_S;
    Double I_ESP_K3xy3z_Px = I_ESP_L4xy3z_S+ABX*I_ESP_K3xy3z_S;
    Double I_ESP_K3x4z_Px = I_ESP_L4x4z_S+ABX*I_ESP_K3x4z_S;
    Double I_ESP_K2x5y_Px = I_ESP_L3x5y_S+ABX*I_ESP_K2x5y_S;
    Double I_ESP_K2x4yz_Px = I_ESP_L3x4yz_S+ABX*I_ESP_K2x4yz_S;
    Double I_ESP_K2x3y2z_Px = I_ESP_L3x3y2z_S+ABX*I_ESP_K2x3y2z_S;
    Double I_ESP_K2x2y3z_Px = I_ESP_L3x2y3z_S+ABX*I_ESP_K2x2y3z_S;
    Double I_ESP_K2xy4z_Px = I_ESP_L3xy4z_S+ABX*I_ESP_K2xy4z_S;
    Double I_ESP_K2x5z_Px = I_ESP_L3x5z_S+ABX*I_ESP_K2x5z_S;
    Double I_ESP_Kx6y_Px = I_ESP_L2x6y_S+ABX*I_ESP_Kx6y_S;
    Double I_ESP_Kx5yz_Px = I_ESP_L2x5yz_S+ABX*I_ESP_Kx5yz_S;
    Double I_ESP_Kx4y2z_Px = I_ESP_L2x4y2z_S+ABX*I_ESP_Kx4y2z_S;
    Double I_ESP_Kx3y3z_Px = I_ESP_L2x3y3z_S+ABX*I_ESP_Kx3y3z_S;
    Double I_ESP_Kx2y4z_Px = I_ESP_L2x2y4z_S+ABX*I_ESP_Kx2y4z_S;
    Double I_ESP_Kxy5z_Px = I_ESP_L2xy5z_S+ABX*I_ESP_Kxy5z_S;
    Double I_ESP_Kx6z_Px = I_ESP_L2x6z_S+ABX*I_ESP_Kx6z_S;
    Double I_ESP_K7y_Px = I_ESP_Lx7y_S+ABX*I_ESP_K7y_S;
    Double I_ESP_K6yz_Px = I_ESP_Lx6yz_S+ABX*I_ESP_K6yz_S;
    Double I_ESP_K5y2z_Px = I_ESP_Lx5y2z_S+ABX*I_ESP_K5y2z_S;
    Double I_ESP_K4y3z_Px = I_ESP_Lx4y3z_S+ABX*I_ESP_K4y3z_S;
    Double I_ESP_K3y4z_Px = I_ESP_Lx3y4z_S+ABX*I_ESP_K3y4z_S;
    Double I_ESP_K2y5z_Px = I_ESP_Lx2y5z_S+ABX*I_ESP_K2y5z_S;
    Double I_ESP_Ky6z_Px = I_ESP_Lxy6z_S+ABX*I_ESP_Ky6z_S;
    Double I_ESP_K7z_Px = I_ESP_Lx7z_S+ABX*I_ESP_K7z_S;
    Double I_ESP_K6xy_Py = I_ESP_L6x2y_S+ABY*I_ESP_K6xy_S;
    Double I_ESP_K6xz_Py = I_ESP_L6xyz_S+ABY*I_ESP_K6xz_S;
    Double I_ESP_K5x2y_Py = I_ESP_L5x3y_S+ABY*I_ESP_K5x2y_S;
    Double I_ESP_K5xyz_Py = I_ESP_L5x2yz_S+ABY*I_ESP_K5xyz_S;
    Double I_ESP_K5x2z_Py = I_ESP_L5xy2z_S+ABY*I_ESP_K5x2z_S;
    Double I_ESP_K4x3y_Py = I_ESP_L4x4y_S+ABY*I_ESP_K4x3y_S;
    Double I_ESP_K4x2yz_Py = I_ESP_L4x3yz_S+ABY*I_ESP_K4x2yz_S;
    Double I_ESP_K4xy2z_Py = I_ESP_L4x2y2z_S+ABY*I_ESP_K4xy2z_S;
    Double I_ESP_K4x3z_Py = I_ESP_L4xy3z_S+ABY*I_ESP_K4x3z_S;
    Double I_ESP_K3x4y_Py = I_ESP_L3x5y_S+ABY*I_ESP_K3x4y_S;
    Double I_ESP_K3x3yz_Py = I_ESP_L3x4yz_S+ABY*I_ESP_K3x3yz_S;
    Double I_ESP_K3x2y2z_Py = I_ESP_L3x3y2z_S+ABY*I_ESP_K3x2y2z_S;
    Double I_ESP_K3xy3z_Py = I_ESP_L3x2y3z_S+ABY*I_ESP_K3xy3z_S;
    Double I_ESP_K3x4z_Py = I_ESP_L3xy4z_S+ABY*I_ESP_K3x4z_S;
    Double I_ESP_K2x5y_Py = I_ESP_L2x6y_S+ABY*I_ESP_K2x5y_S;
    Double I_ESP_K2x4yz_Py = I_ESP_L2x5yz_S+ABY*I_ESP_K2x4yz_S;
    Double I_ESP_K2x3y2z_Py = I_ESP_L2x4y2z_S+ABY*I_ESP_K2x3y2z_S;
    Double I_ESP_K2x2y3z_Py = I_ESP_L2x3y3z_S+ABY*I_ESP_K2x2y3z_S;
    Double I_ESP_K2xy4z_Py = I_ESP_L2x2y4z_S+ABY*I_ESP_K2xy4z_S;
    Double I_ESP_K2x5z_Py = I_ESP_L2xy5z_S+ABY*I_ESP_K2x5z_S;
    Double I_ESP_Kx6y_Py = I_ESP_Lx7y_S+ABY*I_ESP_Kx6y_S;
    Double I_ESP_Kx5yz_Py = I_ESP_Lx6yz_S+ABY*I_ESP_Kx5yz_S;
    Double I_ESP_Kx4y2z_Py = I_ESP_Lx5y2z_S+ABY*I_ESP_Kx4y2z_S;
    Double I_ESP_Kx3y3z_Py = I_ESP_Lx4y3z_S+ABY*I_ESP_Kx3y3z_S;
    Double I_ESP_Kx2y4z_Py = I_ESP_Lx3y4z_S+ABY*I_ESP_Kx2y4z_S;
    Double I_ESP_Kxy5z_Py = I_ESP_Lx2y5z_S+ABY*I_ESP_Kxy5z_S;
    Double I_ESP_Kx6z_Py = I_ESP_Lxy6z_S+ABY*I_ESP_Kx6z_S;
    Double I_ESP_K7y_Py = I_ESP_L8y_S+ABY*I_ESP_K7y_S;
    Double I_ESP_K6yz_Py = I_ESP_L7yz_S+ABY*I_ESP_K6yz_S;
    Double I_ESP_K5y2z_Py = I_ESP_L6y2z_S+ABY*I_ESP_K5y2z_S;
    Double I_ESP_K4y3z_Py = I_ESP_L5y3z_S+ABY*I_ESP_K4y3z_S;
    Double I_ESP_K3y4z_Py = I_ESP_L4y4z_S+ABY*I_ESP_K3y4z_S;
    Double I_ESP_K2y5z_Py = I_ESP_L3y5z_S+ABY*I_ESP_K2y5z_S;
    Double I_ESP_Ky6z_Py = I_ESP_L2y6z_S+ABY*I_ESP_Ky6z_S;
    Double I_ESP_K7z_Py = I_ESP_Ly7z_S+ABY*I_ESP_K7z_S;
    Double I_ESP_K6xy_Pz = I_ESP_L6xyz_S+ABZ*I_ESP_K6xy_S;
    Double I_ESP_K6xz_Pz = I_ESP_L6x2z_S+ABZ*I_ESP_K6xz_S;
    Double I_ESP_K5x2y_Pz = I_ESP_L5x2yz_S+ABZ*I_ESP_K5x2y_S;
    Double I_ESP_K5xyz_Pz = I_ESP_L5xy2z_S+ABZ*I_ESP_K5xyz_S;
    Double I_ESP_K5x2z_Pz = I_ESP_L5x3z_S+ABZ*I_ESP_K5x2z_S;
    Double I_ESP_K4x3y_Pz = I_ESP_L4x3yz_S+ABZ*I_ESP_K4x3y_S;
    Double I_ESP_K4x2yz_Pz = I_ESP_L4x2y2z_S+ABZ*I_ESP_K4x2yz_S;
    Double I_ESP_K4xy2z_Pz = I_ESP_L4xy3z_S+ABZ*I_ESP_K4xy2z_S;
    Double I_ESP_K4x3z_Pz = I_ESP_L4x4z_S+ABZ*I_ESP_K4x3z_S;
    Double I_ESP_K3x4y_Pz = I_ESP_L3x4yz_S+ABZ*I_ESP_K3x4y_S;
    Double I_ESP_K3x3yz_Pz = I_ESP_L3x3y2z_S+ABZ*I_ESP_K3x3yz_S;
    Double I_ESP_K3x2y2z_Pz = I_ESP_L3x2y3z_S+ABZ*I_ESP_K3x2y2z_S;
    Double I_ESP_K3xy3z_Pz = I_ESP_L3xy4z_S+ABZ*I_ESP_K3xy3z_S;
    Double I_ESP_K3x4z_Pz = I_ESP_L3x5z_S+ABZ*I_ESP_K3x4z_S;
    Double I_ESP_K2x5y_Pz = I_ESP_L2x5yz_S+ABZ*I_ESP_K2x5y_S;
    Double I_ESP_K2x4yz_Pz = I_ESP_L2x4y2z_S+ABZ*I_ESP_K2x4yz_S;
    Double I_ESP_K2x3y2z_Pz = I_ESP_L2x3y3z_S+ABZ*I_ESP_K2x3y2z_S;
    Double I_ESP_K2x2y3z_Pz = I_ESP_L2x2y4z_S+ABZ*I_ESP_K2x2y3z_S;
    Double I_ESP_K2xy4z_Pz = I_ESP_L2xy5z_S+ABZ*I_ESP_K2xy4z_S;
    Double I_ESP_K2x5z_Pz = I_ESP_L2x6z_S+ABZ*I_ESP_K2x5z_S;
    Double I_ESP_Kx6y_Pz = I_ESP_Lx6yz_S+ABZ*I_ESP_Kx6y_S;
    Double I_ESP_Kx5yz_Pz = I_ESP_Lx5y2z_S+ABZ*I_ESP_Kx5yz_S;
    Double I_ESP_Kx4y2z_Pz = I_ESP_Lx4y3z_S+ABZ*I_ESP_Kx4y2z_S;
    Double I_ESP_Kx3y3z_Pz = I_ESP_Lx3y4z_S+ABZ*I_ESP_Kx3y3z_S;
    Double I_ESP_Kx2y4z_Pz = I_ESP_Lx2y5z_S+ABZ*I_ESP_Kx2y4z_S;
    Double I_ESP_Kxy5z_Pz = I_ESP_Lxy6z_S+ABZ*I_ESP_Kxy5z_S;
    Double I_ESP_Kx6z_Pz = I_ESP_Lx7z_S+ABZ*I_ESP_Kx6z_S;
    Double I_ESP_K6yz_Pz = I_ESP_L6y2z_S+ABZ*I_ESP_K6yz_S;
    Double I_ESP_K5y2z_Pz = I_ESP_L5y3z_S+ABZ*I_ESP_K5y2z_S;
    Double I_ESP_K4y3z_Pz = I_ESP_L4y4z_S+ABZ*I_ESP_K4y3z_S;
    Double I_ESP_K3y4z_Pz = I_ESP_L3y5z_S+ABZ*I_ESP_K3y4z_S;
    Double I_ESP_K2y5z_Pz = I_ESP_L2y6z_S+ABZ*I_ESP_K2y5z_S;
    Double I_ESP_Ky6z_Pz = I_ESP_Ly7z_S+ABZ*I_ESP_Ky6z_S;
    Double I_ESP_K7z_Pz = I_ESP_L8z_S+ABZ*I_ESP_K7z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_I_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 84 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P
     * RHS shell quartet name: SQ_ESP_I_P
     ************************************************************/
    Double I_ESP_I6x_D2x = I_ESP_K7x_Px+ABX*I_ESP_I6x_Px;
    Double I_ESP_I5xy_D2x = I_ESP_K6xy_Px+ABX*I_ESP_I5xy_Px;
    Double I_ESP_I5xz_D2x = I_ESP_K6xz_Px+ABX*I_ESP_I5xz_Px;
    Double I_ESP_I4x2y_D2x = I_ESP_K5x2y_Px+ABX*I_ESP_I4x2y_Px;
    Double I_ESP_I4xyz_D2x = I_ESP_K5xyz_Px+ABX*I_ESP_I4xyz_Px;
    Double I_ESP_I4x2z_D2x = I_ESP_K5x2z_Px+ABX*I_ESP_I4x2z_Px;
    Double I_ESP_I3x3y_D2x = I_ESP_K4x3y_Px+ABX*I_ESP_I3x3y_Px;
    Double I_ESP_I3x2yz_D2x = I_ESP_K4x2yz_Px+ABX*I_ESP_I3x2yz_Px;
    Double I_ESP_I3xy2z_D2x = I_ESP_K4xy2z_Px+ABX*I_ESP_I3xy2z_Px;
    Double I_ESP_I3x3z_D2x = I_ESP_K4x3z_Px+ABX*I_ESP_I3x3z_Px;
    Double I_ESP_I2x4y_D2x = I_ESP_K3x4y_Px+ABX*I_ESP_I2x4y_Px;
    Double I_ESP_I2x3yz_D2x = I_ESP_K3x3yz_Px+ABX*I_ESP_I2x3yz_Px;
    Double I_ESP_I2x2y2z_D2x = I_ESP_K3x2y2z_Px+ABX*I_ESP_I2x2y2z_Px;
    Double I_ESP_I2xy3z_D2x = I_ESP_K3xy3z_Px+ABX*I_ESP_I2xy3z_Px;
    Double I_ESP_I2x4z_D2x = I_ESP_K3x4z_Px+ABX*I_ESP_I2x4z_Px;
    Double I_ESP_Ix5y_D2x = I_ESP_K2x5y_Px+ABX*I_ESP_Ix5y_Px;
    Double I_ESP_Ix4yz_D2x = I_ESP_K2x4yz_Px+ABX*I_ESP_Ix4yz_Px;
    Double I_ESP_Ix3y2z_D2x = I_ESP_K2x3y2z_Px+ABX*I_ESP_Ix3y2z_Px;
    Double I_ESP_Ix2y3z_D2x = I_ESP_K2x2y3z_Px+ABX*I_ESP_Ix2y3z_Px;
    Double I_ESP_Ixy4z_D2x = I_ESP_K2xy4z_Px+ABX*I_ESP_Ixy4z_Px;
    Double I_ESP_Ix5z_D2x = I_ESP_K2x5z_Px+ABX*I_ESP_Ix5z_Px;
    Double I_ESP_I6y_D2x = I_ESP_Kx6y_Px+ABX*I_ESP_I6y_Px;
    Double I_ESP_I5yz_D2x = I_ESP_Kx5yz_Px+ABX*I_ESP_I5yz_Px;
    Double I_ESP_I4y2z_D2x = I_ESP_Kx4y2z_Px+ABX*I_ESP_I4y2z_Px;
    Double I_ESP_I3y3z_D2x = I_ESP_Kx3y3z_Px+ABX*I_ESP_I3y3z_Px;
    Double I_ESP_I2y4z_D2x = I_ESP_Kx2y4z_Px+ABX*I_ESP_I2y4z_Px;
    Double I_ESP_Iy5z_D2x = I_ESP_Kxy5z_Px+ABX*I_ESP_Iy5z_Px;
    Double I_ESP_I6z_D2x = I_ESP_Kx6z_Px+ABX*I_ESP_I6z_Px;
    Double I_ESP_I6x_D2y = I_ESP_K6xy_Py+ABY*I_ESP_I6x_Py;
    Double I_ESP_I5xy_D2y = I_ESP_K5x2y_Py+ABY*I_ESP_I5xy_Py;
    Double I_ESP_I5xz_D2y = I_ESP_K5xyz_Py+ABY*I_ESP_I5xz_Py;
    Double I_ESP_I4x2y_D2y = I_ESP_K4x3y_Py+ABY*I_ESP_I4x2y_Py;
    Double I_ESP_I4xyz_D2y = I_ESP_K4x2yz_Py+ABY*I_ESP_I4xyz_Py;
    Double I_ESP_I4x2z_D2y = I_ESP_K4xy2z_Py+ABY*I_ESP_I4x2z_Py;
    Double I_ESP_I3x3y_D2y = I_ESP_K3x4y_Py+ABY*I_ESP_I3x3y_Py;
    Double I_ESP_I3x2yz_D2y = I_ESP_K3x3yz_Py+ABY*I_ESP_I3x2yz_Py;
    Double I_ESP_I3xy2z_D2y = I_ESP_K3x2y2z_Py+ABY*I_ESP_I3xy2z_Py;
    Double I_ESP_I3x3z_D2y = I_ESP_K3xy3z_Py+ABY*I_ESP_I3x3z_Py;
    Double I_ESP_I2x4y_D2y = I_ESP_K2x5y_Py+ABY*I_ESP_I2x4y_Py;
    Double I_ESP_I2x3yz_D2y = I_ESP_K2x4yz_Py+ABY*I_ESP_I2x3yz_Py;
    Double I_ESP_I2x2y2z_D2y = I_ESP_K2x3y2z_Py+ABY*I_ESP_I2x2y2z_Py;
    Double I_ESP_I2xy3z_D2y = I_ESP_K2x2y3z_Py+ABY*I_ESP_I2xy3z_Py;
    Double I_ESP_I2x4z_D2y = I_ESP_K2xy4z_Py+ABY*I_ESP_I2x4z_Py;
    Double I_ESP_Ix5y_D2y = I_ESP_Kx6y_Py+ABY*I_ESP_Ix5y_Py;
    Double I_ESP_Ix4yz_D2y = I_ESP_Kx5yz_Py+ABY*I_ESP_Ix4yz_Py;
    Double I_ESP_Ix3y2z_D2y = I_ESP_Kx4y2z_Py+ABY*I_ESP_Ix3y2z_Py;
    Double I_ESP_Ix2y3z_D2y = I_ESP_Kx3y3z_Py+ABY*I_ESP_Ix2y3z_Py;
    Double I_ESP_Ixy4z_D2y = I_ESP_Kx2y4z_Py+ABY*I_ESP_Ixy4z_Py;
    Double I_ESP_Ix5z_D2y = I_ESP_Kxy5z_Py+ABY*I_ESP_Ix5z_Py;
    Double I_ESP_I6y_D2y = I_ESP_K7y_Py+ABY*I_ESP_I6y_Py;
    Double I_ESP_I5yz_D2y = I_ESP_K6yz_Py+ABY*I_ESP_I5yz_Py;
    Double I_ESP_I4y2z_D2y = I_ESP_K5y2z_Py+ABY*I_ESP_I4y2z_Py;
    Double I_ESP_I3y3z_D2y = I_ESP_K4y3z_Py+ABY*I_ESP_I3y3z_Py;
    Double I_ESP_I2y4z_D2y = I_ESP_K3y4z_Py+ABY*I_ESP_I2y4z_Py;
    Double I_ESP_Iy5z_D2y = I_ESP_K2y5z_Py+ABY*I_ESP_Iy5z_Py;
    Double I_ESP_I6z_D2y = I_ESP_Ky6z_Py+ABY*I_ESP_I6z_Py;
    Double I_ESP_I6x_D2z = I_ESP_K6xz_Pz+ABZ*I_ESP_I6x_Pz;
    Double I_ESP_I5xy_D2z = I_ESP_K5xyz_Pz+ABZ*I_ESP_I5xy_Pz;
    Double I_ESP_I5xz_D2z = I_ESP_K5x2z_Pz+ABZ*I_ESP_I5xz_Pz;
    Double I_ESP_I4x2y_D2z = I_ESP_K4x2yz_Pz+ABZ*I_ESP_I4x2y_Pz;
    Double I_ESP_I4xyz_D2z = I_ESP_K4xy2z_Pz+ABZ*I_ESP_I4xyz_Pz;
    Double I_ESP_I4x2z_D2z = I_ESP_K4x3z_Pz+ABZ*I_ESP_I4x2z_Pz;
    Double I_ESP_I3x3y_D2z = I_ESP_K3x3yz_Pz+ABZ*I_ESP_I3x3y_Pz;
    Double I_ESP_I3x2yz_D2z = I_ESP_K3x2y2z_Pz+ABZ*I_ESP_I3x2yz_Pz;
    Double I_ESP_I3xy2z_D2z = I_ESP_K3xy3z_Pz+ABZ*I_ESP_I3xy2z_Pz;
    Double I_ESP_I3x3z_D2z = I_ESP_K3x4z_Pz+ABZ*I_ESP_I3x3z_Pz;
    Double I_ESP_I2x4y_D2z = I_ESP_K2x4yz_Pz+ABZ*I_ESP_I2x4y_Pz;
    Double I_ESP_I2x3yz_D2z = I_ESP_K2x3y2z_Pz+ABZ*I_ESP_I2x3yz_Pz;
    Double I_ESP_I2x2y2z_D2z = I_ESP_K2x2y3z_Pz+ABZ*I_ESP_I2x2y2z_Pz;
    Double I_ESP_I2xy3z_D2z = I_ESP_K2xy4z_Pz+ABZ*I_ESP_I2xy3z_Pz;
    Double I_ESP_I2x4z_D2z = I_ESP_K2x5z_Pz+ABZ*I_ESP_I2x4z_Pz;
    Double I_ESP_Ix5y_D2z = I_ESP_Kx5yz_Pz+ABZ*I_ESP_Ix5y_Pz;
    Double I_ESP_Ix4yz_D2z = I_ESP_Kx4y2z_Pz+ABZ*I_ESP_Ix4yz_Pz;
    Double I_ESP_Ix3y2z_D2z = I_ESP_Kx3y3z_Pz+ABZ*I_ESP_Ix3y2z_Pz;
    Double I_ESP_Ix2y3z_D2z = I_ESP_Kx2y4z_Pz+ABZ*I_ESP_Ix2y3z_Pz;
    Double I_ESP_Ixy4z_D2z = I_ESP_Kxy5z_Pz+ABZ*I_ESP_Ixy4z_Pz;
    Double I_ESP_Ix5z_D2z = I_ESP_Kx6z_Pz+ABZ*I_ESP_Ix5z_Pz;
    Double I_ESP_I6y_D2z = I_ESP_K6yz_Pz+ABZ*I_ESP_I6y_Pz;
    Double I_ESP_I5yz_D2z = I_ESP_K5y2z_Pz+ABZ*I_ESP_I5yz_Pz;
    Double I_ESP_I4y2z_D2z = I_ESP_K4y3z_Pz+ABZ*I_ESP_I4y2z_Pz;
    Double I_ESP_I3y3z_D2z = I_ESP_K3y4z_Pz+ABZ*I_ESP_I3y3z_Pz;
    Double I_ESP_I2y4z_D2z = I_ESP_K2y5z_Pz+ABZ*I_ESP_I2y4z_Pz;
    Double I_ESP_Iy5z_D2z = I_ESP_Ky6z_Pz+ABZ*I_ESP_Iy5z_Pz;
    Double I_ESP_I6z_D2z = I_ESP_K7z_Pz+ABZ*I_ESP_I6z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 42 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D
     * RHS shell quartet name: SQ_ESP_H_D
     ************************************************************/
    Double I_ESP_H5x_F3x = I_ESP_I6x_D2x+ABX*I_ESP_H5x_D2x;
    Double I_ESP_H4xy_F3x = I_ESP_I5xy_D2x+ABX*I_ESP_H4xy_D2x;
    Double I_ESP_H4xz_F3x = I_ESP_I5xz_D2x+ABX*I_ESP_H4xz_D2x;
    Double I_ESP_H3x2y_F3x = I_ESP_I4x2y_D2x+ABX*I_ESP_H3x2y_D2x;
    Double I_ESP_H3xyz_F3x = I_ESP_I4xyz_D2x+ABX*I_ESP_H3xyz_D2x;
    Double I_ESP_H3x2z_F3x = I_ESP_I4x2z_D2x+ABX*I_ESP_H3x2z_D2x;
    Double I_ESP_H2x3y_F3x = I_ESP_I3x3y_D2x+ABX*I_ESP_H2x3y_D2x;
    Double I_ESP_H2x2yz_F3x = I_ESP_I3x2yz_D2x+ABX*I_ESP_H2x2yz_D2x;
    Double I_ESP_H2xy2z_F3x = I_ESP_I3xy2z_D2x+ABX*I_ESP_H2xy2z_D2x;
    Double I_ESP_H2x3z_F3x = I_ESP_I3x3z_D2x+ABX*I_ESP_H2x3z_D2x;
    Double I_ESP_Hx4y_F3x = I_ESP_I2x4y_D2x+ABX*I_ESP_Hx4y_D2x;
    Double I_ESP_Hx3yz_F3x = I_ESP_I2x3yz_D2x+ABX*I_ESP_Hx3yz_D2x;
    Double I_ESP_Hx2y2z_F3x = I_ESP_I2x2y2z_D2x+ABX*I_ESP_Hx2y2z_D2x;
    Double I_ESP_Hxy3z_F3x = I_ESP_I2xy3z_D2x+ABX*I_ESP_Hxy3z_D2x;
    Double I_ESP_Hx4z_F3x = I_ESP_I2x4z_D2x+ABX*I_ESP_Hx4z_D2x;
    Double I_ESP_H5y_F3x = I_ESP_Ix5y_D2x+ABX*I_ESP_H5y_D2x;
    Double I_ESP_H4yz_F3x = I_ESP_Ix4yz_D2x+ABX*I_ESP_H4yz_D2x;
    Double I_ESP_H3y2z_F3x = I_ESP_Ix3y2z_D2x+ABX*I_ESP_H3y2z_D2x;
    Double I_ESP_H2y3z_F3x = I_ESP_Ix2y3z_D2x+ABX*I_ESP_H2y3z_D2x;
    Double I_ESP_Hy4z_F3x = I_ESP_Ixy4z_D2x+ABX*I_ESP_Hy4z_D2x;
    Double I_ESP_H5z_F3x = I_ESP_Ix5z_D2x+ABX*I_ESP_H5z_D2x;
    Double I_ESP_H5x_F2xy = I_ESP_I5xy_D2x+ABY*I_ESP_H5x_D2x;
    Double I_ESP_H4xy_F2xy = I_ESP_I4x2y_D2x+ABY*I_ESP_H4xy_D2x;
    Double I_ESP_H4xz_F2xy = I_ESP_I4xyz_D2x+ABY*I_ESP_H4xz_D2x;
    Double I_ESP_H3x2y_F2xy = I_ESP_I3x3y_D2x+ABY*I_ESP_H3x2y_D2x;
    Double I_ESP_H3xyz_F2xy = I_ESP_I3x2yz_D2x+ABY*I_ESP_H3xyz_D2x;
    Double I_ESP_H3x2z_F2xy = I_ESP_I3xy2z_D2x+ABY*I_ESP_H3x2z_D2x;
    Double I_ESP_H2x3y_F2xy = I_ESP_I2x4y_D2x+ABY*I_ESP_H2x3y_D2x;
    Double I_ESP_H2x2yz_F2xy = I_ESP_I2x3yz_D2x+ABY*I_ESP_H2x2yz_D2x;
    Double I_ESP_H2xy2z_F2xy = I_ESP_I2x2y2z_D2x+ABY*I_ESP_H2xy2z_D2x;
    Double I_ESP_H2x3z_F2xy = I_ESP_I2xy3z_D2x+ABY*I_ESP_H2x3z_D2x;
    Double I_ESP_Hx4y_F2xy = I_ESP_Ix5y_D2x+ABY*I_ESP_Hx4y_D2x;
    Double I_ESP_Hx3yz_F2xy = I_ESP_Ix4yz_D2x+ABY*I_ESP_Hx3yz_D2x;
    Double I_ESP_Hx2y2z_F2xy = I_ESP_Ix3y2z_D2x+ABY*I_ESP_Hx2y2z_D2x;
    Double I_ESP_Hxy3z_F2xy = I_ESP_Ix2y3z_D2x+ABY*I_ESP_Hxy3z_D2x;
    Double I_ESP_Hx4z_F2xy = I_ESP_Ixy4z_D2x+ABY*I_ESP_Hx4z_D2x;
    Double I_ESP_H5y_F2xy = I_ESP_I6y_D2x+ABY*I_ESP_H5y_D2x;
    Double I_ESP_H4yz_F2xy = I_ESP_I5yz_D2x+ABY*I_ESP_H4yz_D2x;
    Double I_ESP_H3y2z_F2xy = I_ESP_I4y2z_D2x+ABY*I_ESP_H3y2z_D2x;
    Double I_ESP_H2y3z_F2xy = I_ESP_I3y3z_D2x+ABY*I_ESP_H2y3z_D2x;
    Double I_ESP_Hy4z_F2xy = I_ESP_I2y4z_D2x+ABY*I_ESP_Hy4z_D2x;
    Double I_ESP_H5z_F2xy = I_ESP_Iy5z_D2x+ABY*I_ESP_H5z_D2x;
    Double I_ESP_H5x_F2xz = I_ESP_I5xz_D2x+ABZ*I_ESP_H5x_D2x;
    Double I_ESP_H4xy_F2xz = I_ESP_I4xyz_D2x+ABZ*I_ESP_H4xy_D2x;
    Double I_ESP_H4xz_F2xz = I_ESP_I4x2z_D2x+ABZ*I_ESP_H4xz_D2x;
    Double I_ESP_H3x2y_F2xz = I_ESP_I3x2yz_D2x+ABZ*I_ESP_H3x2y_D2x;
    Double I_ESP_H3xyz_F2xz = I_ESP_I3xy2z_D2x+ABZ*I_ESP_H3xyz_D2x;
    Double I_ESP_H3x2z_F2xz = I_ESP_I3x3z_D2x+ABZ*I_ESP_H3x2z_D2x;
    Double I_ESP_H2x3y_F2xz = I_ESP_I2x3yz_D2x+ABZ*I_ESP_H2x3y_D2x;
    Double I_ESP_H2x2yz_F2xz = I_ESP_I2x2y2z_D2x+ABZ*I_ESP_H2x2yz_D2x;
    Double I_ESP_H2xy2z_F2xz = I_ESP_I2xy3z_D2x+ABZ*I_ESP_H2xy2z_D2x;
    Double I_ESP_H2x3z_F2xz = I_ESP_I2x4z_D2x+ABZ*I_ESP_H2x3z_D2x;
    Double I_ESP_Hx4y_F2xz = I_ESP_Ix4yz_D2x+ABZ*I_ESP_Hx4y_D2x;
    Double I_ESP_Hx3yz_F2xz = I_ESP_Ix3y2z_D2x+ABZ*I_ESP_Hx3yz_D2x;
    Double I_ESP_Hx2y2z_F2xz = I_ESP_Ix2y3z_D2x+ABZ*I_ESP_Hx2y2z_D2x;
    Double I_ESP_Hxy3z_F2xz = I_ESP_Ixy4z_D2x+ABZ*I_ESP_Hxy3z_D2x;
    Double I_ESP_Hx4z_F2xz = I_ESP_Ix5z_D2x+ABZ*I_ESP_Hx4z_D2x;
    Double I_ESP_H5y_F2xz = I_ESP_I5yz_D2x+ABZ*I_ESP_H5y_D2x;
    Double I_ESP_H4yz_F2xz = I_ESP_I4y2z_D2x+ABZ*I_ESP_H4yz_D2x;
    Double I_ESP_H3y2z_F2xz = I_ESP_I3y3z_D2x+ABZ*I_ESP_H3y2z_D2x;
    Double I_ESP_H2y3z_F2xz = I_ESP_I2y4z_D2x+ABZ*I_ESP_H2y3z_D2x;
    Double I_ESP_Hy4z_F2xz = I_ESP_Iy5z_D2x+ABZ*I_ESP_Hy4z_D2x;
    Double I_ESP_H5z_F2xz = I_ESP_I6z_D2x+ABZ*I_ESP_H5z_D2x;
    Double I_ESP_H5x_Fx2y = I_ESP_I6x_D2y+ABX*I_ESP_H5x_D2y;
    Double I_ESP_H4xy_Fx2y = I_ESP_I5xy_D2y+ABX*I_ESP_H4xy_D2y;
    Double I_ESP_H4xz_Fx2y = I_ESP_I5xz_D2y+ABX*I_ESP_H4xz_D2y;
    Double I_ESP_H3x2y_Fx2y = I_ESP_I4x2y_D2y+ABX*I_ESP_H3x2y_D2y;
    Double I_ESP_H3xyz_Fx2y = I_ESP_I4xyz_D2y+ABX*I_ESP_H3xyz_D2y;
    Double I_ESP_H3x2z_Fx2y = I_ESP_I4x2z_D2y+ABX*I_ESP_H3x2z_D2y;
    Double I_ESP_H2x3y_Fx2y = I_ESP_I3x3y_D2y+ABX*I_ESP_H2x3y_D2y;
    Double I_ESP_H2x2yz_Fx2y = I_ESP_I3x2yz_D2y+ABX*I_ESP_H2x2yz_D2y;
    Double I_ESP_H2xy2z_Fx2y = I_ESP_I3xy2z_D2y+ABX*I_ESP_H2xy2z_D2y;
    Double I_ESP_H2x3z_Fx2y = I_ESP_I3x3z_D2y+ABX*I_ESP_H2x3z_D2y;
    Double I_ESP_Hx4y_Fx2y = I_ESP_I2x4y_D2y+ABX*I_ESP_Hx4y_D2y;
    Double I_ESP_Hx3yz_Fx2y = I_ESP_I2x3yz_D2y+ABX*I_ESP_Hx3yz_D2y;
    Double I_ESP_Hx2y2z_Fx2y = I_ESP_I2x2y2z_D2y+ABX*I_ESP_Hx2y2z_D2y;
    Double I_ESP_Hxy3z_Fx2y = I_ESP_I2xy3z_D2y+ABX*I_ESP_Hxy3z_D2y;
    Double I_ESP_Hx4z_Fx2y = I_ESP_I2x4z_D2y+ABX*I_ESP_Hx4z_D2y;
    Double I_ESP_H5y_Fx2y = I_ESP_Ix5y_D2y+ABX*I_ESP_H5y_D2y;
    Double I_ESP_H4yz_Fx2y = I_ESP_Ix4yz_D2y+ABX*I_ESP_H4yz_D2y;
    Double I_ESP_H3y2z_Fx2y = I_ESP_Ix3y2z_D2y+ABX*I_ESP_H3y2z_D2y;
    Double I_ESP_H2y3z_Fx2y = I_ESP_Ix2y3z_D2y+ABX*I_ESP_H2y3z_D2y;
    Double I_ESP_Hy4z_Fx2y = I_ESP_Ixy4z_D2y+ABX*I_ESP_Hy4z_D2y;
    Double I_ESP_H5z_Fx2y = I_ESP_Ix5z_D2y+ABX*I_ESP_H5z_D2y;
    Double I_ESP_H5x_Fx2z = I_ESP_I6x_D2z+ABX*I_ESP_H5x_D2z;
    Double I_ESP_H4xy_Fx2z = I_ESP_I5xy_D2z+ABX*I_ESP_H4xy_D2z;
    Double I_ESP_H4xz_Fx2z = I_ESP_I5xz_D2z+ABX*I_ESP_H4xz_D2z;
    Double I_ESP_H3x2y_Fx2z = I_ESP_I4x2y_D2z+ABX*I_ESP_H3x2y_D2z;
    Double I_ESP_H3xyz_Fx2z = I_ESP_I4xyz_D2z+ABX*I_ESP_H3xyz_D2z;
    Double I_ESP_H3x2z_Fx2z = I_ESP_I4x2z_D2z+ABX*I_ESP_H3x2z_D2z;
    Double I_ESP_H2x3y_Fx2z = I_ESP_I3x3y_D2z+ABX*I_ESP_H2x3y_D2z;
    Double I_ESP_H2x2yz_Fx2z = I_ESP_I3x2yz_D2z+ABX*I_ESP_H2x2yz_D2z;
    Double I_ESP_H2xy2z_Fx2z = I_ESP_I3xy2z_D2z+ABX*I_ESP_H2xy2z_D2z;
    Double I_ESP_H2x3z_Fx2z = I_ESP_I3x3z_D2z+ABX*I_ESP_H2x3z_D2z;
    Double I_ESP_Hx4y_Fx2z = I_ESP_I2x4y_D2z+ABX*I_ESP_Hx4y_D2z;
    Double I_ESP_Hx3yz_Fx2z = I_ESP_I2x3yz_D2z+ABX*I_ESP_Hx3yz_D2z;
    Double I_ESP_Hx2y2z_Fx2z = I_ESP_I2x2y2z_D2z+ABX*I_ESP_Hx2y2z_D2z;
    Double I_ESP_Hxy3z_Fx2z = I_ESP_I2xy3z_D2z+ABX*I_ESP_Hxy3z_D2z;
    Double I_ESP_Hx4z_Fx2z = I_ESP_I2x4z_D2z+ABX*I_ESP_Hx4z_D2z;
    Double I_ESP_H5y_Fx2z = I_ESP_Ix5y_D2z+ABX*I_ESP_H5y_D2z;
    Double I_ESP_H4yz_Fx2z = I_ESP_Ix4yz_D2z+ABX*I_ESP_H4yz_D2z;
    Double I_ESP_H3y2z_Fx2z = I_ESP_Ix3y2z_D2z+ABX*I_ESP_H3y2z_D2z;
    Double I_ESP_H2y3z_Fx2z = I_ESP_Ix2y3z_D2z+ABX*I_ESP_H2y3z_D2z;
    Double I_ESP_Hy4z_Fx2z = I_ESP_Ixy4z_D2z+ABX*I_ESP_Hy4z_D2z;
    Double I_ESP_H5z_Fx2z = I_ESP_Ix5z_D2z+ABX*I_ESP_H5z_D2z;
    Double I_ESP_H5x_F3y = I_ESP_I5xy_D2y+ABY*I_ESP_H5x_D2y;
    Double I_ESP_H4xy_F3y = I_ESP_I4x2y_D2y+ABY*I_ESP_H4xy_D2y;
    Double I_ESP_H4xz_F3y = I_ESP_I4xyz_D2y+ABY*I_ESP_H4xz_D2y;
    Double I_ESP_H3x2y_F3y = I_ESP_I3x3y_D2y+ABY*I_ESP_H3x2y_D2y;
    Double I_ESP_H3xyz_F3y = I_ESP_I3x2yz_D2y+ABY*I_ESP_H3xyz_D2y;
    Double I_ESP_H3x2z_F3y = I_ESP_I3xy2z_D2y+ABY*I_ESP_H3x2z_D2y;
    Double I_ESP_H2x3y_F3y = I_ESP_I2x4y_D2y+ABY*I_ESP_H2x3y_D2y;
    Double I_ESP_H2x2yz_F3y = I_ESP_I2x3yz_D2y+ABY*I_ESP_H2x2yz_D2y;
    Double I_ESP_H2xy2z_F3y = I_ESP_I2x2y2z_D2y+ABY*I_ESP_H2xy2z_D2y;
    Double I_ESP_H2x3z_F3y = I_ESP_I2xy3z_D2y+ABY*I_ESP_H2x3z_D2y;
    Double I_ESP_Hx4y_F3y = I_ESP_Ix5y_D2y+ABY*I_ESP_Hx4y_D2y;
    Double I_ESP_Hx3yz_F3y = I_ESP_Ix4yz_D2y+ABY*I_ESP_Hx3yz_D2y;
    Double I_ESP_Hx2y2z_F3y = I_ESP_Ix3y2z_D2y+ABY*I_ESP_Hx2y2z_D2y;
    Double I_ESP_Hxy3z_F3y = I_ESP_Ix2y3z_D2y+ABY*I_ESP_Hxy3z_D2y;
    Double I_ESP_Hx4z_F3y = I_ESP_Ixy4z_D2y+ABY*I_ESP_Hx4z_D2y;
    Double I_ESP_H5y_F3y = I_ESP_I6y_D2y+ABY*I_ESP_H5y_D2y;
    Double I_ESP_H4yz_F3y = I_ESP_I5yz_D2y+ABY*I_ESP_H4yz_D2y;
    Double I_ESP_H3y2z_F3y = I_ESP_I4y2z_D2y+ABY*I_ESP_H3y2z_D2y;
    Double I_ESP_H2y3z_F3y = I_ESP_I3y3z_D2y+ABY*I_ESP_H2y3z_D2y;
    Double I_ESP_Hy4z_F3y = I_ESP_I2y4z_D2y+ABY*I_ESP_Hy4z_D2y;
    Double I_ESP_H5z_F3y = I_ESP_Iy5z_D2y+ABY*I_ESP_H5z_D2y;
    Double I_ESP_H5x_F2yz = I_ESP_I5xz_D2y+ABZ*I_ESP_H5x_D2y;
    Double I_ESP_H4xy_F2yz = I_ESP_I4xyz_D2y+ABZ*I_ESP_H4xy_D2y;
    Double I_ESP_H4xz_F2yz = I_ESP_I4x2z_D2y+ABZ*I_ESP_H4xz_D2y;
    Double I_ESP_H3x2y_F2yz = I_ESP_I3x2yz_D2y+ABZ*I_ESP_H3x2y_D2y;
    Double I_ESP_H3xyz_F2yz = I_ESP_I3xy2z_D2y+ABZ*I_ESP_H3xyz_D2y;
    Double I_ESP_H3x2z_F2yz = I_ESP_I3x3z_D2y+ABZ*I_ESP_H3x2z_D2y;
    Double I_ESP_H2x3y_F2yz = I_ESP_I2x3yz_D2y+ABZ*I_ESP_H2x3y_D2y;
    Double I_ESP_H2x2yz_F2yz = I_ESP_I2x2y2z_D2y+ABZ*I_ESP_H2x2yz_D2y;
    Double I_ESP_H2xy2z_F2yz = I_ESP_I2xy3z_D2y+ABZ*I_ESP_H2xy2z_D2y;
    Double I_ESP_H2x3z_F2yz = I_ESP_I2x4z_D2y+ABZ*I_ESP_H2x3z_D2y;
    Double I_ESP_Hx4y_F2yz = I_ESP_Ix4yz_D2y+ABZ*I_ESP_Hx4y_D2y;
    Double I_ESP_Hx3yz_F2yz = I_ESP_Ix3y2z_D2y+ABZ*I_ESP_Hx3yz_D2y;
    Double I_ESP_Hx2y2z_F2yz = I_ESP_Ix2y3z_D2y+ABZ*I_ESP_Hx2y2z_D2y;
    Double I_ESP_Hxy3z_F2yz = I_ESP_Ixy4z_D2y+ABZ*I_ESP_Hxy3z_D2y;
    Double I_ESP_Hx4z_F2yz = I_ESP_Ix5z_D2y+ABZ*I_ESP_Hx4z_D2y;
    Double I_ESP_H5y_F2yz = I_ESP_I5yz_D2y+ABZ*I_ESP_H5y_D2y;
    Double I_ESP_H4yz_F2yz = I_ESP_I4y2z_D2y+ABZ*I_ESP_H4yz_D2y;
    Double I_ESP_H3y2z_F2yz = I_ESP_I3y3z_D2y+ABZ*I_ESP_H3y2z_D2y;
    Double I_ESP_H2y3z_F2yz = I_ESP_I2y4z_D2y+ABZ*I_ESP_H2y3z_D2y;
    Double I_ESP_Hy4z_F2yz = I_ESP_Iy5z_D2y+ABZ*I_ESP_Hy4z_D2y;
    Double I_ESP_H5z_F2yz = I_ESP_I6z_D2y+ABZ*I_ESP_H5z_D2y;
    Double I_ESP_H5x_F3z = I_ESP_I5xz_D2z+ABZ*I_ESP_H5x_D2z;
    Double I_ESP_H4xy_F3z = I_ESP_I4xyz_D2z+ABZ*I_ESP_H4xy_D2z;
    Double I_ESP_H4xz_F3z = I_ESP_I4x2z_D2z+ABZ*I_ESP_H4xz_D2z;
    Double I_ESP_H3x2y_F3z = I_ESP_I3x2yz_D2z+ABZ*I_ESP_H3x2y_D2z;
    Double I_ESP_H3xyz_F3z = I_ESP_I3xy2z_D2z+ABZ*I_ESP_H3xyz_D2z;
    Double I_ESP_H3x2z_F3z = I_ESP_I3x3z_D2z+ABZ*I_ESP_H3x2z_D2z;
    Double I_ESP_H2x3y_F3z = I_ESP_I2x3yz_D2z+ABZ*I_ESP_H2x3y_D2z;
    Double I_ESP_H2x2yz_F3z = I_ESP_I2x2y2z_D2z+ABZ*I_ESP_H2x2yz_D2z;
    Double I_ESP_H2xy2z_F3z = I_ESP_I2xy3z_D2z+ABZ*I_ESP_H2xy2z_D2z;
    Double I_ESP_H2x3z_F3z = I_ESP_I2x4z_D2z+ABZ*I_ESP_H2x3z_D2z;
    Double I_ESP_Hx4y_F3z = I_ESP_Ix4yz_D2z+ABZ*I_ESP_Hx4y_D2z;
    Double I_ESP_Hx3yz_F3z = I_ESP_Ix3y2z_D2z+ABZ*I_ESP_Hx3yz_D2z;
    Double I_ESP_Hx2y2z_F3z = I_ESP_Ix2y3z_D2z+ABZ*I_ESP_Hx2y2z_D2z;
    Double I_ESP_Hxy3z_F3z = I_ESP_Ixy4z_D2z+ABZ*I_ESP_Hxy3z_D2z;
    Double I_ESP_Hx4z_F3z = I_ESP_Ix5z_D2z+ABZ*I_ESP_Hx4z_D2z;
    Double I_ESP_H5y_F3z = I_ESP_I5yz_D2z+ABZ*I_ESP_H5y_D2z;
    Double I_ESP_H4yz_F3z = I_ESP_I4y2z_D2z+ABZ*I_ESP_H4yz_D2z;
    Double I_ESP_H3y2z_F3z = I_ESP_I3y3z_D2z+ABZ*I_ESP_H3y2z_D2z;
    Double I_ESP_H2y3z_F3z = I_ESP_I2y4z_D2z+ABZ*I_ESP_H2y3z_D2z;
    Double I_ESP_Hy4z_F3z = I_ESP_Iy5z_D2z+ABZ*I_ESP_Hy4z_D2z;
    Double I_ESP_H5z_F3z = I_ESP_I6z_D2z+ABZ*I_ESP_H5z_D2z;

    /************************************************************
     * shell quartet name: SQ_ESP_L_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_M_S
     * RHS shell quartet name: SQ_ESP_L_S
     ************************************************************/
    Double I_ESP_L8x_Px = I_ESP_M9x_S+ABX*I_ESP_L8x_S;
    Double I_ESP_L7xy_Px = I_ESP_M8xy_S+ABX*I_ESP_L7xy_S;
    Double I_ESP_L7xz_Px = I_ESP_M8xz_S+ABX*I_ESP_L7xz_S;
    Double I_ESP_L6x2y_Px = I_ESP_M7x2y_S+ABX*I_ESP_L6x2y_S;
    Double I_ESP_L6xyz_Px = I_ESP_M7xyz_S+ABX*I_ESP_L6xyz_S;
    Double I_ESP_L6x2z_Px = I_ESP_M7x2z_S+ABX*I_ESP_L6x2z_S;
    Double I_ESP_L5x3y_Px = I_ESP_M6x3y_S+ABX*I_ESP_L5x3y_S;
    Double I_ESP_L5x2yz_Px = I_ESP_M6x2yz_S+ABX*I_ESP_L5x2yz_S;
    Double I_ESP_L5xy2z_Px = I_ESP_M6xy2z_S+ABX*I_ESP_L5xy2z_S;
    Double I_ESP_L5x3z_Px = I_ESP_M6x3z_S+ABX*I_ESP_L5x3z_S;
    Double I_ESP_L4x4y_Px = I_ESP_M5x4y_S+ABX*I_ESP_L4x4y_S;
    Double I_ESP_L4x3yz_Px = I_ESP_M5x3yz_S+ABX*I_ESP_L4x3yz_S;
    Double I_ESP_L4x2y2z_Px = I_ESP_M5x2y2z_S+ABX*I_ESP_L4x2y2z_S;
    Double I_ESP_L4xy3z_Px = I_ESP_M5xy3z_S+ABX*I_ESP_L4xy3z_S;
    Double I_ESP_L4x4z_Px = I_ESP_M5x4z_S+ABX*I_ESP_L4x4z_S;
    Double I_ESP_L3x5y_Px = I_ESP_M4x5y_S+ABX*I_ESP_L3x5y_S;
    Double I_ESP_L3x4yz_Px = I_ESP_M4x4yz_S+ABX*I_ESP_L3x4yz_S;
    Double I_ESP_L3x3y2z_Px = I_ESP_M4x3y2z_S+ABX*I_ESP_L3x3y2z_S;
    Double I_ESP_L3x2y3z_Px = I_ESP_M4x2y3z_S+ABX*I_ESP_L3x2y3z_S;
    Double I_ESP_L3xy4z_Px = I_ESP_M4xy4z_S+ABX*I_ESP_L3xy4z_S;
    Double I_ESP_L3x5z_Px = I_ESP_M4x5z_S+ABX*I_ESP_L3x5z_S;
    Double I_ESP_L2x6y_Px = I_ESP_M3x6y_S+ABX*I_ESP_L2x6y_S;
    Double I_ESP_L2x5yz_Px = I_ESP_M3x5yz_S+ABX*I_ESP_L2x5yz_S;
    Double I_ESP_L2x4y2z_Px = I_ESP_M3x4y2z_S+ABX*I_ESP_L2x4y2z_S;
    Double I_ESP_L2x3y3z_Px = I_ESP_M3x3y3z_S+ABX*I_ESP_L2x3y3z_S;
    Double I_ESP_L2x2y4z_Px = I_ESP_M3x2y4z_S+ABX*I_ESP_L2x2y4z_S;
    Double I_ESP_L2xy5z_Px = I_ESP_M3xy5z_S+ABX*I_ESP_L2xy5z_S;
    Double I_ESP_L2x6z_Px = I_ESP_M3x6z_S+ABX*I_ESP_L2x6z_S;
    Double I_ESP_Lx7y_Px = I_ESP_M2x7y_S+ABX*I_ESP_Lx7y_S;
    Double I_ESP_Lx6yz_Px = I_ESP_M2x6yz_S+ABX*I_ESP_Lx6yz_S;
    Double I_ESP_Lx5y2z_Px = I_ESP_M2x5y2z_S+ABX*I_ESP_Lx5y2z_S;
    Double I_ESP_Lx4y3z_Px = I_ESP_M2x4y3z_S+ABX*I_ESP_Lx4y3z_S;
    Double I_ESP_Lx3y4z_Px = I_ESP_M2x3y4z_S+ABX*I_ESP_Lx3y4z_S;
    Double I_ESP_Lx2y5z_Px = I_ESP_M2x2y5z_S+ABX*I_ESP_Lx2y5z_S;
    Double I_ESP_Lxy6z_Px = I_ESP_M2xy6z_S+ABX*I_ESP_Lxy6z_S;
    Double I_ESP_Lx7z_Px = I_ESP_M2x7z_S+ABX*I_ESP_Lx7z_S;
    Double I_ESP_L6x2y_Py = I_ESP_M6x3y_S+ABY*I_ESP_L6x2y_S;
    Double I_ESP_L6xyz_Py = I_ESP_M6x2yz_S+ABY*I_ESP_L6xyz_S;
    Double I_ESP_L5x3y_Py = I_ESP_M5x4y_S+ABY*I_ESP_L5x3y_S;
    Double I_ESP_L5x2yz_Py = I_ESP_M5x3yz_S+ABY*I_ESP_L5x2yz_S;
    Double I_ESP_L5xy2z_Py = I_ESP_M5x2y2z_S+ABY*I_ESP_L5xy2z_S;
    Double I_ESP_L4x4y_Py = I_ESP_M4x5y_S+ABY*I_ESP_L4x4y_S;
    Double I_ESP_L4x3yz_Py = I_ESP_M4x4yz_S+ABY*I_ESP_L4x3yz_S;
    Double I_ESP_L4x2y2z_Py = I_ESP_M4x3y2z_S+ABY*I_ESP_L4x2y2z_S;
    Double I_ESP_L4xy3z_Py = I_ESP_M4x2y3z_S+ABY*I_ESP_L4xy3z_S;
    Double I_ESP_L3x5y_Py = I_ESP_M3x6y_S+ABY*I_ESP_L3x5y_S;
    Double I_ESP_L3x4yz_Py = I_ESP_M3x5yz_S+ABY*I_ESP_L3x4yz_S;
    Double I_ESP_L3x3y2z_Py = I_ESP_M3x4y2z_S+ABY*I_ESP_L3x3y2z_S;
    Double I_ESP_L3x2y3z_Py = I_ESP_M3x3y3z_S+ABY*I_ESP_L3x2y3z_S;
    Double I_ESP_L3xy4z_Py = I_ESP_M3x2y4z_S+ABY*I_ESP_L3xy4z_S;
    Double I_ESP_L2x6y_Py = I_ESP_M2x7y_S+ABY*I_ESP_L2x6y_S;
    Double I_ESP_L2x5yz_Py = I_ESP_M2x6yz_S+ABY*I_ESP_L2x5yz_S;
    Double I_ESP_L2x4y2z_Py = I_ESP_M2x5y2z_S+ABY*I_ESP_L2x4y2z_S;
    Double I_ESP_L2x3y3z_Py = I_ESP_M2x4y3z_S+ABY*I_ESP_L2x3y3z_S;
    Double I_ESP_L2x2y4z_Py = I_ESP_M2x3y4z_S+ABY*I_ESP_L2x2y4z_S;
    Double I_ESP_L2xy5z_Py = I_ESP_M2x2y5z_S+ABY*I_ESP_L2xy5z_S;
    Double I_ESP_Lx7y_Py = I_ESP_Mx8y_S+ABY*I_ESP_Lx7y_S;
    Double I_ESP_Lx6yz_Py = I_ESP_Mx7yz_S+ABY*I_ESP_Lx6yz_S;
    Double I_ESP_Lx5y2z_Py = I_ESP_Mx6y2z_S+ABY*I_ESP_Lx5y2z_S;
    Double I_ESP_Lx4y3z_Py = I_ESP_Mx5y3z_S+ABY*I_ESP_Lx4y3z_S;
    Double I_ESP_Lx3y4z_Py = I_ESP_Mx4y4z_S+ABY*I_ESP_Lx3y4z_S;
    Double I_ESP_Lx2y5z_Py = I_ESP_Mx3y5z_S+ABY*I_ESP_Lx2y5z_S;
    Double I_ESP_Lxy6z_Py = I_ESP_Mx2y6z_S+ABY*I_ESP_Lxy6z_S;
    Double I_ESP_L8y_Py = I_ESP_M9y_S+ABY*I_ESP_L8y_S;
    Double I_ESP_L7yz_Py = I_ESP_M8yz_S+ABY*I_ESP_L7yz_S;
    Double I_ESP_L6y2z_Py = I_ESP_M7y2z_S+ABY*I_ESP_L6y2z_S;
    Double I_ESP_L5y3z_Py = I_ESP_M6y3z_S+ABY*I_ESP_L5y3z_S;
    Double I_ESP_L4y4z_Py = I_ESP_M5y4z_S+ABY*I_ESP_L4y4z_S;
    Double I_ESP_L3y5z_Py = I_ESP_M4y5z_S+ABY*I_ESP_L3y5z_S;
    Double I_ESP_L2y6z_Py = I_ESP_M3y6z_S+ABY*I_ESP_L2y6z_S;
    Double I_ESP_Ly7z_Py = I_ESP_M2y7z_S+ABY*I_ESP_Ly7z_S;
    Double I_ESP_L6xyz_Pz = I_ESP_M6xy2z_S+ABZ*I_ESP_L6xyz_S;
    Double I_ESP_L6x2z_Pz = I_ESP_M6x3z_S+ABZ*I_ESP_L6x2z_S;
    Double I_ESP_L5x2yz_Pz = I_ESP_M5x2y2z_S+ABZ*I_ESP_L5x2yz_S;
    Double I_ESP_L5xy2z_Pz = I_ESP_M5xy3z_S+ABZ*I_ESP_L5xy2z_S;
    Double I_ESP_L5x3z_Pz = I_ESP_M5x4z_S+ABZ*I_ESP_L5x3z_S;
    Double I_ESP_L4x3yz_Pz = I_ESP_M4x3y2z_S+ABZ*I_ESP_L4x3yz_S;
    Double I_ESP_L4x2y2z_Pz = I_ESP_M4x2y3z_S+ABZ*I_ESP_L4x2y2z_S;
    Double I_ESP_L4xy3z_Pz = I_ESP_M4xy4z_S+ABZ*I_ESP_L4xy3z_S;
    Double I_ESP_L4x4z_Pz = I_ESP_M4x5z_S+ABZ*I_ESP_L4x4z_S;
    Double I_ESP_L3x4yz_Pz = I_ESP_M3x4y2z_S+ABZ*I_ESP_L3x4yz_S;
    Double I_ESP_L3x3y2z_Pz = I_ESP_M3x3y3z_S+ABZ*I_ESP_L3x3y2z_S;
    Double I_ESP_L3x2y3z_Pz = I_ESP_M3x2y4z_S+ABZ*I_ESP_L3x2y3z_S;
    Double I_ESP_L3xy4z_Pz = I_ESP_M3xy5z_S+ABZ*I_ESP_L3xy4z_S;
    Double I_ESP_L3x5z_Pz = I_ESP_M3x6z_S+ABZ*I_ESP_L3x5z_S;
    Double I_ESP_L2x5yz_Pz = I_ESP_M2x5y2z_S+ABZ*I_ESP_L2x5yz_S;
    Double I_ESP_L2x4y2z_Pz = I_ESP_M2x4y3z_S+ABZ*I_ESP_L2x4y2z_S;
    Double I_ESP_L2x3y3z_Pz = I_ESP_M2x3y4z_S+ABZ*I_ESP_L2x3y3z_S;
    Double I_ESP_L2x2y4z_Pz = I_ESP_M2x2y5z_S+ABZ*I_ESP_L2x2y4z_S;
    Double I_ESP_L2xy5z_Pz = I_ESP_M2xy6z_S+ABZ*I_ESP_L2xy5z_S;
    Double I_ESP_L2x6z_Pz = I_ESP_M2x7z_S+ABZ*I_ESP_L2x6z_S;
    Double I_ESP_Lx6yz_Pz = I_ESP_Mx6y2z_S+ABZ*I_ESP_Lx6yz_S;
    Double I_ESP_Lx5y2z_Pz = I_ESP_Mx5y3z_S+ABZ*I_ESP_Lx5y2z_S;
    Double I_ESP_Lx4y3z_Pz = I_ESP_Mx4y4z_S+ABZ*I_ESP_Lx4y3z_S;
    Double I_ESP_Lx3y4z_Pz = I_ESP_Mx3y5z_S+ABZ*I_ESP_Lx3y4z_S;
    Double I_ESP_Lx2y5z_Pz = I_ESP_Mx2y6z_S+ABZ*I_ESP_Lx2y5z_S;
    Double I_ESP_Lxy6z_Pz = I_ESP_Mxy7z_S+ABZ*I_ESP_Lxy6z_S;
    Double I_ESP_Lx7z_Pz = I_ESP_Mx8z_S+ABZ*I_ESP_Lx7z_S;
    Double I_ESP_L6y2z_Pz = I_ESP_M6y3z_S+ABZ*I_ESP_L6y2z_S;
    Double I_ESP_L5y3z_Pz = I_ESP_M5y4z_S+ABZ*I_ESP_L5y3z_S;
    Double I_ESP_L4y4z_Pz = I_ESP_M4y5z_S+ABZ*I_ESP_L4y4z_S;
    Double I_ESP_L3y5z_Pz = I_ESP_M3y6z_S+ABZ*I_ESP_L3y5z_S;
    Double I_ESP_L2y6z_Pz = I_ESP_M2y7z_S+ABZ*I_ESP_L2y6z_S;
    Double I_ESP_Ly7z_Pz = I_ESP_My8z_S+ABZ*I_ESP_Ly7z_S;
    Double I_ESP_L8z_Pz = I_ESP_M9z_S+ABZ*I_ESP_L8z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_K_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 111 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_P
     * RHS shell quartet name: SQ_ESP_K_P
     ************************************************************/
    Double I_ESP_K7x_D2x = I_ESP_L8x_Px+ABX*I_ESP_K7x_Px;
    Double I_ESP_K6xy_D2x = I_ESP_L7xy_Px+ABX*I_ESP_K6xy_Px;
    Double I_ESP_K6xz_D2x = I_ESP_L7xz_Px+ABX*I_ESP_K6xz_Px;
    Double I_ESP_K5x2y_D2x = I_ESP_L6x2y_Px+ABX*I_ESP_K5x2y_Px;
    Double I_ESP_K5xyz_D2x = I_ESP_L6xyz_Px+ABX*I_ESP_K5xyz_Px;
    Double I_ESP_K5x2z_D2x = I_ESP_L6x2z_Px+ABX*I_ESP_K5x2z_Px;
    Double I_ESP_K4x3y_D2x = I_ESP_L5x3y_Px+ABX*I_ESP_K4x3y_Px;
    Double I_ESP_K4x2yz_D2x = I_ESP_L5x2yz_Px+ABX*I_ESP_K4x2yz_Px;
    Double I_ESP_K4xy2z_D2x = I_ESP_L5xy2z_Px+ABX*I_ESP_K4xy2z_Px;
    Double I_ESP_K4x3z_D2x = I_ESP_L5x3z_Px+ABX*I_ESP_K4x3z_Px;
    Double I_ESP_K3x4y_D2x = I_ESP_L4x4y_Px+ABX*I_ESP_K3x4y_Px;
    Double I_ESP_K3x3yz_D2x = I_ESP_L4x3yz_Px+ABX*I_ESP_K3x3yz_Px;
    Double I_ESP_K3x2y2z_D2x = I_ESP_L4x2y2z_Px+ABX*I_ESP_K3x2y2z_Px;
    Double I_ESP_K3xy3z_D2x = I_ESP_L4xy3z_Px+ABX*I_ESP_K3xy3z_Px;
    Double I_ESP_K3x4z_D2x = I_ESP_L4x4z_Px+ABX*I_ESP_K3x4z_Px;
    Double I_ESP_K2x5y_D2x = I_ESP_L3x5y_Px+ABX*I_ESP_K2x5y_Px;
    Double I_ESP_K2x4yz_D2x = I_ESP_L3x4yz_Px+ABX*I_ESP_K2x4yz_Px;
    Double I_ESP_K2x3y2z_D2x = I_ESP_L3x3y2z_Px+ABX*I_ESP_K2x3y2z_Px;
    Double I_ESP_K2x2y3z_D2x = I_ESP_L3x2y3z_Px+ABX*I_ESP_K2x2y3z_Px;
    Double I_ESP_K2xy4z_D2x = I_ESP_L3xy4z_Px+ABX*I_ESP_K2xy4z_Px;
    Double I_ESP_K2x5z_D2x = I_ESP_L3x5z_Px+ABX*I_ESP_K2x5z_Px;
    Double I_ESP_Kx6y_D2x = I_ESP_L2x6y_Px+ABX*I_ESP_Kx6y_Px;
    Double I_ESP_Kx5yz_D2x = I_ESP_L2x5yz_Px+ABX*I_ESP_Kx5yz_Px;
    Double I_ESP_Kx4y2z_D2x = I_ESP_L2x4y2z_Px+ABX*I_ESP_Kx4y2z_Px;
    Double I_ESP_Kx3y3z_D2x = I_ESP_L2x3y3z_Px+ABX*I_ESP_Kx3y3z_Px;
    Double I_ESP_Kx2y4z_D2x = I_ESP_L2x2y4z_Px+ABX*I_ESP_Kx2y4z_Px;
    Double I_ESP_Kxy5z_D2x = I_ESP_L2xy5z_Px+ABX*I_ESP_Kxy5z_Px;
    Double I_ESP_Kx6z_D2x = I_ESP_L2x6z_Px+ABX*I_ESP_Kx6z_Px;
    Double I_ESP_K7y_D2x = I_ESP_Lx7y_Px+ABX*I_ESP_K7y_Px;
    Double I_ESP_K6yz_D2x = I_ESP_Lx6yz_Px+ABX*I_ESP_K6yz_Px;
    Double I_ESP_K5y2z_D2x = I_ESP_Lx5y2z_Px+ABX*I_ESP_K5y2z_Px;
    Double I_ESP_K4y3z_D2x = I_ESP_Lx4y3z_Px+ABX*I_ESP_K4y3z_Px;
    Double I_ESP_K3y4z_D2x = I_ESP_Lx3y4z_Px+ABX*I_ESP_K3y4z_Px;
    Double I_ESP_K2y5z_D2x = I_ESP_Lx2y5z_Px+ABX*I_ESP_K2y5z_Px;
    Double I_ESP_Ky6z_D2x = I_ESP_Lxy6z_Px+ABX*I_ESP_Ky6z_Px;
    Double I_ESP_K7z_D2x = I_ESP_Lx7z_Px+ABX*I_ESP_K7z_Px;
    Double I_ESP_K6xy_D2y = I_ESP_L6x2y_Py+ABY*I_ESP_K6xy_Py;
    Double I_ESP_K6xz_D2y = I_ESP_L6xyz_Py+ABY*I_ESP_K6xz_Py;
    Double I_ESP_K5x2y_D2y = I_ESP_L5x3y_Py+ABY*I_ESP_K5x2y_Py;
    Double I_ESP_K5xyz_D2y = I_ESP_L5x2yz_Py+ABY*I_ESP_K5xyz_Py;
    Double I_ESP_K5x2z_D2y = I_ESP_L5xy2z_Py+ABY*I_ESP_K5x2z_Py;
    Double I_ESP_K4x3y_D2y = I_ESP_L4x4y_Py+ABY*I_ESP_K4x3y_Py;
    Double I_ESP_K4x2yz_D2y = I_ESP_L4x3yz_Py+ABY*I_ESP_K4x2yz_Py;
    Double I_ESP_K4xy2z_D2y = I_ESP_L4x2y2z_Py+ABY*I_ESP_K4xy2z_Py;
    Double I_ESP_K4x3z_D2y = I_ESP_L4xy3z_Py+ABY*I_ESP_K4x3z_Py;
    Double I_ESP_K3x4y_D2y = I_ESP_L3x5y_Py+ABY*I_ESP_K3x4y_Py;
    Double I_ESP_K3x3yz_D2y = I_ESP_L3x4yz_Py+ABY*I_ESP_K3x3yz_Py;
    Double I_ESP_K3x2y2z_D2y = I_ESP_L3x3y2z_Py+ABY*I_ESP_K3x2y2z_Py;
    Double I_ESP_K3xy3z_D2y = I_ESP_L3x2y3z_Py+ABY*I_ESP_K3xy3z_Py;
    Double I_ESP_K3x4z_D2y = I_ESP_L3xy4z_Py+ABY*I_ESP_K3x4z_Py;
    Double I_ESP_K2x5y_D2y = I_ESP_L2x6y_Py+ABY*I_ESP_K2x5y_Py;
    Double I_ESP_K2x4yz_D2y = I_ESP_L2x5yz_Py+ABY*I_ESP_K2x4yz_Py;
    Double I_ESP_K2x3y2z_D2y = I_ESP_L2x4y2z_Py+ABY*I_ESP_K2x3y2z_Py;
    Double I_ESP_K2x2y3z_D2y = I_ESP_L2x3y3z_Py+ABY*I_ESP_K2x2y3z_Py;
    Double I_ESP_K2xy4z_D2y = I_ESP_L2x2y4z_Py+ABY*I_ESP_K2xy4z_Py;
    Double I_ESP_K2x5z_D2y = I_ESP_L2xy5z_Py+ABY*I_ESP_K2x5z_Py;
    Double I_ESP_Kx6y_D2y = I_ESP_Lx7y_Py+ABY*I_ESP_Kx6y_Py;
    Double I_ESP_Kx5yz_D2y = I_ESP_Lx6yz_Py+ABY*I_ESP_Kx5yz_Py;
    Double I_ESP_Kx4y2z_D2y = I_ESP_Lx5y2z_Py+ABY*I_ESP_Kx4y2z_Py;
    Double I_ESP_Kx3y3z_D2y = I_ESP_Lx4y3z_Py+ABY*I_ESP_Kx3y3z_Py;
    Double I_ESP_Kx2y4z_D2y = I_ESP_Lx3y4z_Py+ABY*I_ESP_Kx2y4z_Py;
    Double I_ESP_Kxy5z_D2y = I_ESP_Lx2y5z_Py+ABY*I_ESP_Kxy5z_Py;
    Double I_ESP_Kx6z_D2y = I_ESP_Lxy6z_Py+ABY*I_ESP_Kx6z_Py;
    Double I_ESP_K7y_D2y = I_ESP_L8y_Py+ABY*I_ESP_K7y_Py;
    Double I_ESP_K6yz_D2y = I_ESP_L7yz_Py+ABY*I_ESP_K6yz_Py;
    Double I_ESP_K5y2z_D2y = I_ESP_L6y2z_Py+ABY*I_ESP_K5y2z_Py;
    Double I_ESP_K4y3z_D2y = I_ESP_L5y3z_Py+ABY*I_ESP_K4y3z_Py;
    Double I_ESP_K3y4z_D2y = I_ESP_L4y4z_Py+ABY*I_ESP_K3y4z_Py;
    Double I_ESP_K2y5z_D2y = I_ESP_L3y5z_Py+ABY*I_ESP_K2y5z_Py;
    Double I_ESP_Ky6z_D2y = I_ESP_L2y6z_Py+ABY*I_ESP_Ky6z_Py;
    Double I_ESP_K7z_D2y = I_ESP_Ly7z_Py+ABY*I_ESP_K7z_Py;
    Double I_ESP_K6xy_D2z = I_ESP_L6xyz_Pz+ABZ*I_ESP_K6xy_Pz;
    Double I_ESP_K6xz_D2z = I_ESP_L6x2z_Pz+ABZ*I_ESP_K6xz_Pz;
    Double I_ESP_K5x2y_D2z = I_ESP_L5x2yz_Pz+ABZ*I_ESP_K5x2y_Pz;
    Double I_ESP_K5xyz_D2z = I_ESP_L5xy2z_Pz+ABZ*I_ESP_K5xyz_Pz;
    Double I_ESP_K5x2z_D2z = I_ESP_L5x3z_Pz+ABZ*I_ESP_K5x2z_Pz;
    Double I_ESP_K4x3y_D2z = I_ESP_L4x3yz_Pz+ABZ*I_ESP_K4x3y_Pz;
    Double I_ESP_K4x2yz_D2z = I_ESP_L4x2y2z_Pz+ABZ*I_ESP_K4x2yz_Pz;
    Double I_ESP_K4xy2z_D2z = I_ESP_L4xy3z_Pz+ABZ*I_ESP_K4xy2z_Pz;
    Double I_ESP_K4x3z_D2z = I_ESP_L4x4z_Pz+ABZ*I_ESP_K4x3z_Pz;
    Double I_ESP_K3x4y_D2z = I_ESP_L3x4yz_Pz+ABZ*I_ESP_K3x4y_Pz;
    Double I_ESP_K3x3yz_D2z = I_ESP_L3x3y2z_Pz+ABZ*I_ESP_K3x3yz_Pz;
    Double I_ESP_K3x2y2z_D2z = I_ESP_L3x2y3z_Pz+ABZ*I_ESP_K3x2y2z_Pz;
    Double I_ESP_K3xy3z_D2z = I_ESP_L3xy4z_Pz+ABZ*I_ESP_K3xy3z_Pz;
    Double I_ESP_K3x4z_D2z = I_ESP_L3x5z_Pz+ABZ*I_ESP_K3x4z_Pz;
    Double I_ESP_K2x5y_D2z = I_ESP_L2x5yz_Pz+ABZ*I_ESP_K2x5y_Pz;
    Double I_ESP_K2x4yz_D2z = I_ESP_L2x4y2z_Pz+ABZ*I_ESP_K2x4yz_Pz;
    Double I_ESP_K2x3y2z_D2z = I_ESP_L2x3y3z_Pz+ABZ*I_ESP_K2x3y2z_Pz;
    Double I_ESP_K2x2y3z_D2z = I_ESP_L2x2y4z_Pz+ABZ*I_ESP_K2x2y3z_Pz;
    Double I_ESP_K2xy4z_D2z = I_ESP_L2xy5z_Pz+ABZ*I_ESP_K2xy4z_Pz;
    Double I_ESP_K2x5z_D2z = I_ESP_L2x6z_Pz+ABZ*I_ESP_K2x5z_Pz;
    Double I_ESP_Kx6y_D2z = I_ESP_Lx6yz_Pz+ABZ*I_ESP_Kx6y_Pz;
    Double I_ESP_Kx5yz_D2z = I_ESP_Lx5y2z_Pz+ABZ*I_ESP_Kx5yz_Pz;
    Double I_ESP_Kx4y2z_D2z = I_ESP_Lx4y3z_Pz+ABZ*I_ESP_Kx4y2z_Pz;
    Double I_ESP_Kx3y3z_D2z = I_ESP_Lx3y4z_Pz+ABZ*I_ESP_Kx3y3z_Pz;
    Double I_ESP_Kx2y4z_D2z = I_ESP_Lx2y5z_Pz+ABZ*I_ESP_Kx2y4z_Pz;
    Double I_ESP_Kxy5z_D2z = I_ESP_Lxy6z_Pz+ABZ*I_ESP_Kxy5z_Pz;
    Double I_ESP_Kx6z_D2z = I_ESP_Lx7z_Pz+ABZ*I_ESP_Kx6z_Pz;
    Double I_ESP_K6yz_D2z = I_ESP_L6y2z_Pz+ABZ*I_ESP_K6yz_Pz;
    Double I_ESP_K5y2z_D2z = I_ESP_L5y3z_Pz+ABZ*I_ESP_K5y2z_Pz;
    Double I_ESP_K4y3z_D2z = I_ESP_L4y4z_Pz+ABZ*I_ESP_K4y3z_Pz;
    Double I_ESP_K3y4z_D2z = I_ESP_L3y5z_Pz+ABZ*I_ESP_K3y4z_Pz;
    Double I_ESP_K2y5z_D2z = I_ESP_L2y6z_Pz+ABZ*I_ESP_K2y5z_Pz;
    Double I_ESP_Ky6z_D2z = I_ESP_Ly7z_Pz+ABZ*I_ESP_Ky6z_Pz;
    Double I_ESP_K7z_D2z = I_ESP_L8z_Pz+ABZ*I_ESP_K7z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_I_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 85 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_D
     * RHS shell quartet name: SQ_ESP_I_D
     ************************************************************/
    Double I_ESP_I6x_F3x = I_ESP_K7x_D2x+ABX*I_ESP_I6x_D2x;
    Double I_ESP_I5xy_F3x = I_ESP_K6xy_D2x+ABX*I_ESP_I5xy_D2x;
    Double I_ESP_I5xz_F3x = I_ESP_K6xz_D2x+ABX*I_ESP_I5xz_D2x;
    Double I_ESP_I4x2y_F3x = I_ESP_K5x2y_D2x+ABX*I_ESP_I4x2y_D2x;
    Double I_ESP_I4xyz_F3x = I_ESP_K5xyz_D2x+ABX*I_ESP_I4xyz_D2x;
    Double I_ESP_I4x2z_F3x = I_ESP_K5x2z_D2x+ABX*I_ESP_I4x2z_D2x;
    Double I_ESP_I3x3y_F3x = I_ESP_K4x3y_D2x+ABX*I_ESP_I3x3y_D2x;
    Double I_ESP_I3x2yz_F3x = I_ESP_K4x2yz_D2x+ABX*I_ESP_I3x2yz_D2x;
    Double I_ESP_I3xy2z_F3x = I_ESP_K4xy2z_D2x+ABX*I_ESP_I3xy2z_D2x;
    Double I_ESP_I3x3z_F3x = I_ESP_K4x3z_D2x+ABX*I_ESP_I3x3z_D2x;
    Double I_ESP_I2x4y_F3x = I_ESP_K3x4y_D2x+ABX*I_ESP_I2x4y_D2x;
    Double I_ESP_I2x3yz_F3x = I_ESP_K3x3yz_D2x+ABX*I_ESP_I2x3yz_D2x;
    Double I_ESP_I2x2y2z_F3x = I_ESP_K3x2y2z_D2x+ABX*I_ESP_I2x2y2z_D2x;
    Double I_ESP_I2xy3z_F3x = I_ESP_K3xy3z_D2x+ABX*I_ESP_I2xy3z_D2x;
    Double I_ESP_I2x4z_F3x = I_ESP_K3x4z_D2x+ABX*I_ESP_I2x4z_D2x;
    Double I_ESP_Ix5y_F3x = I_ESP_K2x5y_D2x+ABX*I_ESP_Ix5y_D2x;
    Double I_ESP_Ix4yz_F3x = I_ESP_K2x4yz_D2x+ABX*I_ESP_Ix4yz_D2x;
    Double I_ESP_Ix3y2z_F3x = I_ESP_K2x3y2z_D2x+ABX*I_ESP_Ix3y2z_D2x;
    Double I_ESP_Ix2y3z_F3x = I_ESP_K2x2y3z_D2x+ABX*I_ESP_Ix2y3z_D2x;
    Double I_ESP_Ixy4z_F3x = I_ESP_K2xy4z_D2x+ABX*I_ESP_Ixy4z_D2x;
    Double I_ESP_Ix5z_F3x = I_ESP_K2x5z_D2x+ABX*I_ESP_Ix5z_D2x;
    Double I_ESP_I6y_F3x = I_ESP_Kx6y_D2x+ABX*I_ESP_I6y_D2x;
    Double I_ESP_I5yz_F3x = I_ESP_Kx5yz_D2x+ABX*I_ESP_I5yz_D2x;
    Double I_ESP_I4y2z_F3x = I_ESP_Kx4y2z_D2x+ABX*I_ESP_I4y2z_D2x;
    Double I_ESP_I3y3z_F3x = I_ESP_Kx3y3z_D2x+ABX*I_ESP_I3y3z_D2x;
    Double I_ESP_I2y4z_F3x = I_ESP_Kx2y4z_D2x+ABX*I_ESP_I2y4z_D2x;
    Double I_ESP_Iy5z_F3x = I_ESP_Kxy5z_D2x+ABX*I_ESP_Iy5z_D2x;
    Double I_ESP_I6z_F3x = I_ESP_Kx6z_D2x+ABX*I_ESP_I6z_D2x;
    Double I_ESP_I5xy_F2xy = I_ESP_K5x2y_D2x+ABY*I_ESP_I5xy_D2x;
    Double I_ESP_I5xz_F2xy = I_ESP_K5xyz_D2x+ABY*I_ESP_I5xz_D2x;
    Double I_ESP_I4x2y_F2xy = I_ESP_K4x3y_D2x+ABY*I_ESP_I4x2y_D2x;
    Double I_ESP_I4xyz_F2xy = I_ESP_K4x2yz_D2x+ABY*I_ESP_I4xyz_D2x;
    Double I_ESP_I4x2z_F2xy = I_ESP_K4xy2z_D2x+ABY*I_ESP_I4x2z_D2x;
    Double I_ESP_I3x3y_F2xy = I_ESP_K3x4y_D2x+ABY*I_ESP_I3x3y_D2x;
    Double I_ESP_I3x2yz_F2xy = I_ESP_K3x3yz_D2x+ABY*I_ESP_I3x2yz_D2x;
    Double I_ESP_I3xy2z_F2xy = I_ESP_K3x2y2z_D2x+ABY*I_ESP_I3xy2z_D2x;
    Double I_ESP_I3x3z_F2xy = I_ESP_K3xy3z_D2x+ABY*I_ESP_I3x3z_D2x;
    Double I_ESP_I2x4y_F2xy = I_ESP_K2x5y_D2x+ABY*I_ESP_I2x4y_D2x;
    Double I_ESP_I2x3yz_F2xy = I_ESP_K2x4yz_D2x+ABY*I_ESP_I2x3yz_D2x;
    Double I_ESP_I2x2y2z_F2xy = I_ESP_K2x3y2z_D2x+ABY*I_ESP_I2x2y2z_D2x;
    Double I_ESP_I2xy3z_F2xy = I_ESP_K2x2y3z_D2x+ABY*I_ESP_I2xy3z_D2x;
    Double I_ESP_I2x4z_F2xy = I_ESP_K2xy4z_D2x+ABY*I_ESP_I2x4z_D2x;
    Double I_ESP_Ix5y_F2xy = I_ESP_Kx6y_D2x+ABY*I_ESP_Ix5y_D2x;
    Double I_ESP_Ix4yz_F2xy = I_ESP_Kx5yz_D2x+ABY*I_ESP_Ix4yz_D2x;
    Double I_ESP_Ix3y2z_F2xy = I_ESP_Kx4y2z_D2x+ABY*I_ESP_Ix3y2z_D2x;
    Double I_ESP_Ix2y3z_F2xy = I_ESP_Kx3y3z_D2x+ABY*I_ESP_Ix2y3z_D2x;
    Double I_ESP_Ixy4z_F2xy = I_ESP_Kx2y4z_D2x+ABY*I_ESP_Ixy4z_D2x;
    Double I_ESP_Ix5z_F2xy = I_ESP_Kxy5z_D2x+ABY*I_ESP_Ix5z_D2x;
    Double I_ESP_I6y_F2xy = I_ESP_K7y_D2x+ABY*I_ESP_I6y_D2x;
    Double I_ESP_I5yz_F2xy = I_ESP_K6yz_D2x+ABY*I_ESP_I5yz_D2x;
    Double I_ESP_I4y2z_F2xy = I_ESP_K5y2z_D2x+ABY*I_ESP_I4y2z_D2x;
    Double I_ESP_I3y3z_F2xy = I_ESP_K4y3z_D2x+ABY*I_ESP_I3y3z_D2x;
    Double I_ESP_I2y4z_F2xy = I_ESP_K3y4z_D2x+ABY*I_ESP_I2y4z_D2x;
    Double I_ESP_Iy5z_F2xy = I_ESP_K2y5z_D2x+ABY*I_ESP_Iy5z_D2x;
    Double I_ESP_I6z_F2xy = I_ESP_Ky6z_D2x+ABY*I_ESP_I6z_D2x;
    Double I_ESP_I5xz_F2xz = I_ESP_K5x2z_D2x+ABZ*I_ESP_I5xz_D2x;
    Double I_ESP_I4xyz_F2xz = I_ESP_K4xy2z_D2x+ABZ*I_ESP_I4xyz_D2x;
    Double I_ESP_I4x2z_F2xz = I_ESP_K4x3z_D2x+ABZ*I_ESP_I4x2z_D2x;
    Double I_ESP_I3x2yz_F2xz = I_ESP_K3x2y2z_D2x+ABZ*I_ESP_I3x2yz_D2x;
    Double I_ESP_I3xy2z_F2xz = I_ESP_K3xy3z_D2x+ABZ*I_ESP_I3xy2z_D2x;
    Double I_ESP_I3x3z_F2xz = I_ESP_K3x4z_D2x+ABZ*I_ESP_I3x3z_D2x;
    Double I_ESP_I2x3yz_F2xz = I_ESP_K2x3y2z_D2x+ABZ*I_ESP_I2x3yz_D2x;
    Double I_ESP_I2x2y2z_F2xz = I_ESP_K2x2y3z_D2x+ABZ*I_ESP_I2x2y2z_D2x;
    Double I_ESP_I2xy3z_F2xz = I_ESP_K2xy4z_D2x+ABZ*I_ESP_I2xy3z_D2x;
    Double I_ESP_I2x4z_F2xz = I_ESP_K2x5z_D2x+ABZ*I_ESP_I2x4z_D2x;
    Double I_ESP_Ix4yz_F2xz = I_ESP_Kx4y2z_D2x+ABZ*I_ESP_Ix4yz_D2x;
    Double I_ESP_Ix3y2z_F2xz = I_ESP_Kx3y3z_D2x+ABZ*I_ESP_Ix3y2z_D2x;
    Double I_ESP_Ix2y3z_F2xz = I_ESP_Kx2y4z_D2x+ABZ*I_ESP_Ix2y3z_D2x;
    Double I_ESP_Ixy4z_F2xz = I_ESP_Kxy5z_D2x+ABZ*I_ESP_Ixy4z_D2x;
    Double I_ESP_Ix5z_F2xz = I_ESP_Kx6z_D2x+ABZ*I_ESP_Ix5z_D2x;
    Double I_ESP_I5yz_F2xz = I_ESP_K5y2z_D2x+ABZ*I_ESP_I5yz_D2x;
    Double I_ESP_I4y2z_F2xz = I_ESP_K4y3z_D2x+ABZ*I_ESP_I4y2z_D2x;
    Double I_ESP_I3y3z_F2xz = I_ESP_K3y4z_D2x+ABZ*I_ESP_I3y3z_D2x;
    Double I_ESP_I2y4z_F2xz = I_ESP_K2y5z_D2x+ABZ*I_ESP_I2y4z_D2x;
    Double I_ESP_Iy5z_F2xz = I_ESP_Ky6z_D2x+ABZ*I_ESP_Iy5z_D2x;
    Double I_ESP_I6z_F2xz = I_ESP_K7z_D2x+ABZ*I_ESP_I6z_D2x;
    Double I_ESP_I5xz_Fx2y = I_ESP_K6xz_D2y+ABX*I_ESP_I5xz_D2y;
    Double I_ESP_I4xyz_Fx2y = I_ESP_K5xyz_D2y+ABX*I_ESP_I4xyz_D2y;
    Double I_ESP_I4x2z_Fx2y = I_ESP_K5x2z_D2y+ABX*I_ESP_I4x2z_D2y;
    Double I_ESP_I3x2yz_Fx2y = I_ESP_K4x2yz_D2y+ABX*I_ESP_I3x2yz_D2y;
    Double I_ESP_I3xy2z_Fx2y = I_ESP_K4xy2z_D2y+ABX*I_ESP_I3xy2z_D2y;
    Double I_ESP_I3x3z_Fx2y = I_ESP_K4x3z_D2y+ABX*I_ESP_I3x3z_D2y;
    Double I_ESP_I2x3yz_Fx2y = I_ESP_K3x3yz_D2y+ABX*I_ESP_I2x3yz_D2y;
    Double I_ESP_I2x2y2z_Fx2y = I_ESP_K3x2y2z_D2y+ABX*I_ESP_I2x2y2z_D2y;
    Double I_ESP_I2xy3z_Fx2y = I_ESP_K3xy3z_D2y+ABX*I_ESP_I2xy3z_D2y;
    Double I_ESP_I2x4z_Fx2y = I_ESP_K3x4z_D2y+ABX*I_ESP_I2x4z_D2y;
    Double I_ESP_Ix4yz_Fx2y = I_ESP_K2x4yz_D2y+ABX*I_ESP_Ix4yz_D2y;
    Double I_ESP_Ix3y2z_Fx2y = I_ESP_K2x3y2z_D2y+ABX*I_ESP_Ix3y2z_D2y;
    Double I_ESP_Ix2y3z_Fx2y = I_ESP_K2x2y3z_D2y+ABX*I_ESP_Ix2y3z_D2y;
    Double I_ESP_Ixy4z_Fx2y = I_ESP_K2xy4z_D2y+ABX*I_ESP_Ixy4z_D2y;
    Double I_ESP_Ix5z_Fx2y = I_ESP_K2x5z_D2y+ABX*I_ESP_Ix5z_D2y;
    Double I_ESP_I5yz_Fx2y = I_ESP_Kx5yz_D2y+ABX*I_ESP_I5yz_D2y;
    Double I_ESP_I4y2z_Fx2y = I_ESP_Kx4y2z_D2y+ABX*I_ESP_I4y2z_D2y;
    Double I_ESP_I3y3z_Fx2y = I_ESP_Kx3y3z_D2y+ABX*I_ESP_I3y3z_D2y;
    Double I_ESP_I2y4z_Fx2y = I_ESP_Kx2y4z_D2y+ABX*I_ESP_I2y4z_D2y;
    Double I_ESP_Iy5z_Fx2y = I_ESP_Kxy5z_D2y+ABX*I_ESP_Iy5z_D2y;
    Double I_ESP_I6z_Fx2y = I_ESP_Kx6z_D2y+ABX*I_ESP_I6z_D2y;
    Double I_ESP_I5xy_Fx2z = I_ESP_K6xy_D2z+ABX*I_ESP_I5xy_D2z;
    Double I_ESP_I4x2y_Fx2z = I_ESP_K5x2y_D2z+ABX*I_ESP_I4x2y_D2z;
    Double I_ESP_I4xyz_Fx2z = I_ESP_K5xyz_D2z+ABX*I_ESP_I4xyz_D2z;
    Double I_ESP_I3x3y_Fx2z = I_ESP_K4x3y_D2z+ABX*I_ESP_I3x3y_D2z;
    Double I_ESP_I3x2yz_Fx2z = I_ESP_K4x2yz_D2z+ABX*I_ESP_I3x2yz_D2z;
    Double I_ESP_I3xy2z_Fx2z = I_ESP_K4xy2z_D2z+ABX*I_ESP_I3xy2z_D2z;
    Double I_ESP_I2x4y_Fx2z = I_ESP_K3x4y_D2z+ABX*I_ESP_I2x4y_D2z;
    Double I_ESP_I2x3yz_Fx2z = I_ESP_K3x3yz_D2z+ABX*I_ESP_I2x3yz_D2z;
    Double I_ESP_I2x2y2z_Fx2z = I_ESP_K3x2y2z_D2z+ABX*I_ESP_I2x2y2z_D2z;
    Double I_ESP_I2xy3z_Fx2z = I_ESP_K3xy3z_D2z+ABX*I_ESP_I2xy3z_D2z;
    Double I_ESP_Ix5y_Fx2z = I_ESP_K2x5y_D2z+ABX*I_ESP_Ix5y_D2z;
    Double I_ESP_Ix4yz_Fx2z = I_ESP_K2x4yz_D2z+ABX*I_ESP_Ix4yz_D2z;
    Double I_ESP_Ix3y2z_Fx2z = I_ESP_K2x3y2z_D2z+ABX*I_ESP_Ix3y2z_D2z;
    Double I_ESP_Ix2y3z_Fx2z = I_ESP_K2x2y3z_D2z+ABX*I_ESP_Ix2y3z_D2z;
    Double I_ESP_Ixy4z_Fx2z = I_ESP_K2xy4z_D2z+ABX*I_ESP_Ixy4z_D2z;
    Double I_ESP_I6y_Fx2z = I_ESP_Kx6y_D2z+ABX*I_ESP_I6y_D2z;
    Double I_ESP_I5yz_Fx2z = I_ESP_Kx5yz_D2z+ABX*I_ESP_I5yz_D2z;
    Double I_ESP_I4y2z_Fx2z = I_ESP_Kx4y2z_D2z+ABX*I_ESP_I4y2z_D2z;
    Double I_ESP_I3y3z_Fx2z = I_ESP_Kx3y3z_D2z+ABX*I_ESP_I3y3z_D2z;
    Double I_ESP_I2y4z_Fx2z = I_ESP_Kx2y4z_D2z+ABX*I_ESP_I2y4z_D2z;
    Double I_ESP_Iy5z_Fx2z = I_ESP_Kxy5z_D2z+ABX*I_ESP_Iy5z_D2z;
    Double I_ESP_I6x_F3y = I_ESP_K6xy_D2y+ABY*I_ESP_I6x_D2y;
    Double I_ESP_I5xy_F3y = I_ESP_K5x2y_D2y+ABY*I_ESP_I5xy_D2y;
    Double I_ESP_I5xz_F3y = I_ESP_K5xyz_D2y+ABY*I_ESP_I5xz_D2y;
    Double I_ESP_I4x2y_F3y = I_ESP_K4x3y_D2y+ABY*I_ESP_I4x2y_D2y;
    Double I_ESP_I4xyz_F3y = I_ESP_K4x2yz_D2y+ABY*I_ESP_I4xyz_D2y;
    Double I_ESP_I4x2z_F3y = I_ESP_K4xy2z_D2y+ABY*I_ESP_I4x2z_D2y;
    Double I_ESP_I3x3y_F3y = I_ESP_K3x4y_D2y+ABY*I_ESP_I3x3y_D2y;
    Double I_ESP_I3x2yz_F3y = I_ESP_K3x3yz_D2y+ABY*I_ESP_I3x2yz_D2y;
    Double I_ESP_I3xy2z_F3y = I_ESP_K3x2y2z_D2y+ABY*I_ESP_I3xy2z_D2y;
    Double I_ESP_I3x3z_F3y = I_ESP_K3xy3z_D2y+ABY*I_ESP_I3x3z_D2y;
    Double I_ESP_I2x4y_F3y = I_ESP_K2x5y_D2y+ABY*I_ESP_I2x4y_D2y;
    Double I_ESP_I2x3yz_F3y = I_ESP_K2x4yz_D2y+ABY*I_ESP_I2x3yz_D2y;
    Double I_ESP_I2x2y2z_F3y = I_ESP_K2x3y2z_D2y+ABY*I_ESP_I2x2y2z_D2y;
    Double I_ESP_I2xy3z_F3y = I_ESP_K2x2y3z_D2y+ABY*I_ESP_I2xy3z_D2y;
    Double I_ESP_I2x4z_F3y = I_ESP_K2xy4z_D2y+ABY*I_ESP_I2x4z_D2y;
    Double I_ESP_Ix5y_F3y = I_ESP_Kx6y_D2y+ABY*I_ESP_Ix5y_D2y;
    Double I_ESP_Ix4yz_F3y = I_ESP_Kx5yz_D2y+ABY*I_ESP_Ix4yz_D2y;
    Double I_ESP_Ix3y2z_F3y = I_ESP_Kx4y2z_D2y+ABY*I_ESP_Ix3y2z_D2y;
    Double I_ESP_Ix2y3z_F3y = I_ESP_Kx3y3z_D2y+ABY*I_ESP_Ix2y3z_D2y;
    Double I_ESP_Ixy4z_F3y = I_ESP_Kx2y4z_D2y+ABY*I_ESP_Ixy4z_D2y;
    Double I_ESP_Ix5z_F3y = I_ESP_Kxy5z_D2y+ABY*I_ESP_Ix5z_D2y;
    Double I_ESP_I6y_F3y = I_ESP_K7y_D2y+ABY*I_ESP_I6y_D2y;
    Double I_ESP_I5yz_F3y = I_ESP_K6yz_D2y+ABY*I_ESP_I5yz_D2y;
    Double I_ESP_I4y2z_F3y = I_ESP_K5y2z_D2y+ABY*I_ESP_I4y2z_D2y;
    Double I_ESP_I3y3z_F3y = I_ESP_K4y3z_D2y+ABY*I_ESP_I3y3z_D2y;
    Double I_ESP_I2y4z_F3y = I_ESP_K3y4z_D2y+ABY*I_ESP_I2y4z_D2y;
    Double I_ESP_Iy5z_F3y = I_ESP_K2y5z_D2y+ABY*I_ESP_Iy5z_D2y;
    Double I_ESP_I6z_F3y = I_ESP_Ky6z_D2y+ABY*I_ESP_I6z_D2y;
    Double I_ESP_I5xz_F2yz = I_ESP_K5x2z_D2y+ABZ*I_ESP_I5xz_D2y;
    Double I_ESP_I4xyz_F2yz = I_ESP_K4xy2z_D2y+ABZ*I_ESP_I4xyz_D2y;
    Double I_ESP_I4x2z_F2yz = I_ESP_K4x3z_D2y+ABZ*I_ESP_I4x2z_D2y;
    Double I_ESP_I3x2yz_F2yz = I_ESP_K3x2y2z_D2y+ABZ*I_ESP_I3x2yz_D2y;
    Double I_ESP_I3xy2z_F2yz = I_ESP_K3xy3z_D2y+ABZ*I_ESP_I3xy2z_D2y;
    Double I_ESP_I3x3z_F2yz = I_ESP_K3x4z_D2y+ABZ*I_ESP_I3x3z_D2y;
    Double I_ESP_I2x3yz_F2yz = I_ESP_K2x3y2z_D2y+ABZ*I_ESP_I2x3yz_D2y;
    Double I_ESP_I2x2y2z_F2yz = I_ESP_K2x2y3z_D2y+ABZ*I_ESP_I2x2y2z_D2y;
    Double I_ESP_I2xy3z_F2yz = I_ESP_K2xy4z_D2y+ABZ*I_ESP_I2xy3z_D2y;
    Double I_ESP_I2x4z_F2yz = I_ESP_K2x5z_D2y+ABZ*I_ESP_I2x4z_D2y;
    Double I_ESP_Ix4yz_F2yz = I_ESP_Kx4y2z_D2y+ABZ*I_ESP_Ix4yz_D2y;
    Double I_ESP_Ix3y2z_F2yz = I_ESP_Kx3y3z_D2y+ABZ*I_ESP_Ix3y2z_D2y;
    Double I_ESP_Ix2y3z_F2yz = I_ESP_Kx2y4z_D2y+ABZ*I_ESP_Ix2y3z_D2y;
    Double I_ESP_Ixy4z_F2yz = I_ESP_Kxy5z_D2y+ABZ*I_ESP_Ixy4z_D2y;
    Double I_ESP_Ix5z_F2yz = I_ESP_Kx6z_D2y+ABZ*I_ESP_Ix5z_D2y;
    Double I_ESP_I5yz_F2yz = I_ESP_K5y2z_D2y+ABZ*I_ESP_I5yz_D2y;
    Double I_ESP_I4y2z_F2yz = I_ESP_K4y3z_D2y+ABZ*I_ESP_I4y2z_D2y;
    Double I_ESP_I3y3z_F2yz = I_ESP_K3y4z_D2y+ABZ*I_ESP_I3y3z_D2y;
    Double I_ESP_I2y4z_F2yz = I_ESP_K2y5z_D2y+ABZ*I_ESP_I2y4z_D2y;
    Double I_ESP_Iy5z_F2yz = I_ESP_Ky6z_D2y+ABZ*I_ESP_Iy5z_D2y;
    Double I_ESP_I6z_F2yz = I_ESP_K7z_D2y+ABZ*I_ESP_I6z_D2y;
    Double I_ESP_I6x_F3z = I_ESP_K6xz_D2z+ABZ*I_ESP_I6x_D2z;
    Double I_ESP_I5xy_F3z = I_ESP_K5xyz_D2z+ABZ*I_ESP_I5xy_D2z;
    Double I_ESP_I5xz_F3z = I_ESP_K5x2z_D2z+ABZ*I_ESP_I5xz_D2z;
    Double I_ESP_I4x2y_F3z = I_ESP_K4x2yz_D2z+ABZ*I_ESP_I4x2y_D2z;
    Double I_ESP_I4xyz_F3z = I_ESP_K4xy2z_D2z+ABZ*I_ESP_I4xyz_D2z;
    Double I_ESP_I4x2z_F3z = I_ESP_K4x3z_D2z+ABZ*I_ESP_I4x2z_D2z;
    Double I_ESP_I3x3y_F3z = I_ESP_K3x3yz_D2z+ABZ*I_ESP_I3x3y_D2z;
    Double I_ESP_I3x2yz_F3z = I_ESP_K3x2y2z_D2z+ABZ*I_ESP_I3x2yz_D2z;
    Double I_ESP_I3xy2z_F3z = I_ESP_K3xy3z_D2z+ABZ*I_ESP_I3xy2z_D2z;
    Double I_ESP_I3x3z_F3z = I_ESP_K3x4z_D2z+ABZ*I_ESP_I3x3z_D2z;
    Double I_ESP_I2x4y_F3z = I_ESP_K2x4yz_D2z+ABZ*I_ESP_I2x4y_D2z;
    Double I_ESP_I2x3yz_F3z = I_ESP_K2x3y2z_D2z+ABZ*I_ESP_I2x3yz_D2z;
    Double I_ESP_I2x2y2z_F3z = I_ESP_K2x2y3z_D2z+ABZ*I_ESP_I2x2y2z_D2z;
    Double I_ESP_I2xy3z_F3z = I_ESP_K2xy4z_D2z+ABZ*I_ESP_I2xy3z_D2z;
    Double I_ESP_I2x4z_F3z = I_ESP_K2x5z_D2z+ABZ*I_ESP_I2x4z_D2z;
    Double I_ESP_Ix5y_F3z = I_ESP_Kx5yz_D2z+ABZ*I_ESP_Ix5y_D2z;
    Double I_ESP_Ix4yz_F3z = I_ESP_Kx4y2z_D2z+ABZ*I_ESP_Ix4yz_D2z;
    Double I_ESP_Ix3y2z_F3z = I_ESP_Kx3y3z_D2z+ABZ*I_ESP_Ix3y2z_D2z;
    Double I_ESP_Ix2y3z_F3z = I_ESP_Kx2y4z_D2z+ABZ*I_ESP_Ix2y3z_D2z;
    Double I_ESP_Ixy4z_F3z = I_ESP_Kxy5z_D2z+ABZ*I_ESP_Ixy4z_D2z;
    Double I_ESP_Ix5z_F3z = I_ESP_Kx6z_D2z+ABZ*I_ESP_Ix5z_D2z;
    Double I_ESP_I6y_F3z = I_ESP_K6yz_D2z+ABZ*I_ESP_I6y_D2z;
    Double I_ESP_I5yz_F3z = I_ESP_K5y2z_D2z+ABZ*I_ESP_I5yz_D2z;
    Double I_ESP_I4y2z_F3z = I_ESP_K4y3z_D2z+ABZ*I_ESP_I4y2z_D2z;
    Double I_ESP_I3y3z_F3z = I_ESP_K3y4z_D2z+ABZ*I_ESP_I3y3z_D2z;
    Double I_ESP_I2y4z_F3z = I_ESP_K2y5z_D2z+ABZ*I_ESP_I2y4z_D2z;
    Double I_ESP_Iy5z_F3z = I_ESP_Ky6z_D2z+ABZ*I_ESP_Iy5z_D2z;
    Double I_ESP_I6z_F3z = I_ESP_K7z_D2z+ABZ*I_ESP_I6z_D2z;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F
     * RHS shell quartet name: SQ_ESP_H_F
     ************************************************************/
    abcd[iGrid*315+0] = I_ESP_I6x_F3x+ABX*I_ESP_H5x_F3x;
    abcd[iGrid*315+1] = I_ESP_I5xy_F3x+ABX*I_ESP_H4xy_F3x;
    abcd[iGrid*315+2] = I_ESP_I5xz_F3x+ABX*I_ESP_H4xz_F3x;
    abcd[iGrid*315+3] = I_ESP_I4x2y_F3x+ABX*I_ESP_H3x2y_F3x;
    abcd[iGrid*315+4] = I_ESP_I4xyz_F3x+ABX*I_ESP_H3xyz_F3x;
    abcd[iGrid*315+5] = I_ESP_I4x2z_F3x+ABX*I_ESP_H3x2z_F3x;
    abcd[iGrid*315+6] = I_ESP_I3x3y_F3x+ABX*I_ESP_H2x3y_F3x;
    abcd[iGrid*315+7] = I_ESP_I3x2yz_F3x+ABX*I_ESP_H2x2yz_F3x;
    abcd[iGrid*315+8] = I_ESP_I3xy2z_F3x+ABX*I_ESP_H2xy2z_F3x;
    abcd[iGrid*315+9] = I_ESP_I3x3z_F3x+ABX*I_ESP_H2x3z_F3x;
    abcd[iGrid*315+10] = I_ESP_I2x4y_F3x+ABX*I_ESP_Hx4y_F3x;
    abcd[iGrid*315+11] = I_ESP_I2x3yz_F3x+ABX*I_ESP_Hx3yz_F3x;
    abcd[iGrid*315+12] = I_ESP_I2x2y2z_F3x+ABX*I_ESP_Hx2y2z_F3x;
    abcd[iGrid*315+13] = I_ESP_I2xy3z_F3x+ABX*I_ESP_Hxy3z_F3x;
    abcd[iGrid*315+14] = I_ESP_I2x4z_F3x+ABX*I_ESP_Hx4z_F3x;
    abcd[iGrid*315+15] = I_ESP_Ix5y_F3x+ABX*I_ESP_H5y_F3x;
    abcd[iGrid*315+16] = I_ESP_Ix4yz_F3x+ABX*I_ESP_H4yz_F3x;
    abcd[iGrid*315+17] = I_ESP_Ix3y2z_F3x+ABX*I_ESP_H3y2z_F3x;
    abcd[iGrid*315+18] = I_ESP_Ix2y3z_F3x+ABX*I_ESP_H2y3z_F3x;
    abcd[iGrid*315+19] = I_ESP_Ixy4z_F3x+ABX*I_ESP_Hy4z_F3x;
    abcd[iGrid*315+20] = I_ESP_Ix5z_F3x+ABX*I_ESP_H5z_F3x;
    abcd[iGrid*315+21] = I_ESP_I5xy_F3x+ABY*I_ESP_H5x_F3x;
    abcd[iGrid*315+22] = I_ESP_I4x2y_F3x+ABY*I_ESP_H4xy_F3x;
    abcd[iGrid*315+23] = I_ESP_I4xyz_F3x+ABY*I_ESP_H4xz_F3x;
    abcd[iGrid*315+24] = I_ESP_I3x3y_F3x+ABY*I_ESP_H3x2y_F3x;
    abcd[iGrid*315+25] = I_ESP_I3x2yz_F3x+ABY*I_ESP_H3xyz_F3x;
    abcd[iGrid*315+26] = I_ESP_I3xy2z_F3x+ABY*I_ESP_H3x2z_F3x;
    abcd[iGrid*315+27] = I_ESP_I2x4y_F3x+ABY*I_ESP_H2x3y_F3x;
    abcd[iGrid*315+28] = I_ESP_I2x3yz_F3x+ABY*I_ESP_H2x2yz_F3x;
    abcd[iGrid*315+29] = I_ESP_I2x2y2z_F3x+ABY*I_ESP_H2xy2z_F3x;
    abcd[iGrid*315+30] = I_ESP_I2xy3z_F3x+ABY*I_ESP_H2x3z_F3x;
    abcd[iGrid*315+31] = I_ESP_Ix5y_F3x+ABY*I_ESP_Hx4y_F3x;
    abcd[iGrid*315+32] = I_ESP_Ix4yz_F3x+ABY*I_ESP_Hx3yz_F3x;
    abcd[iGrid*315+33] = I_ESP_Ix3y2z_F3x+ABY*I_ESP_Hx2y2z_F3x;
    abcd[iGrid*315+34] = I_ESP_Ix2y3z_F3x+ABY*I_ESP_Hxy3z_F3x;
    abcd[iGrid*315+35] = I_ESP_Ixy4z_F3x+ABY*I_ESP_Hx4z_F3x;
    abcd[iGrid*315+36] = I_ESP_I6y_F3x+ABY*I_ESP_H5y_F3x;
    abcd[iGrid*315+37] = I_ESP_I5yz_F3x+ABY*I_ESP_H4yz_F3x;
    abcd[iGrid*315+38] = I_ESP_I4y2z_F3x+ABY*I_ESP_H3y2z_F3x;
    abcd[iGrid*315+39] = I_ESP_I3y3z_F3x+ABY*I_ESP_H2y3z_F3x;
    abcd[iGrid*315+40] = I_ESP_I2y4z_F3x+ABY*I_ESP_Hy4z_F3x;
    abcd[iGrid*315+41] = I_ESP_Iy5z_F3x+ABY*I_ESP_H5z_F3x;
    abcd[iGrid*315+42] = I_ESP_I5xz_F3x+ABZ*I_ESP_H5x_F3x;
    abcd[iGrid*315+43] = I_ESP_I4xyz_F3x+ABZ*I_ESP_H4xy_F3x;
    abcd[iGrid*315+44] = I_ESP_I4x2z_F3x+ABZ*I_ESP_H4xz_F3x;
    abcd[iGrid*315+45] = I_ESP_I3x2yz_F3x+ABZ*I_ESP_H3x2y_F3x;
    abcd[iGrid*315+46] = I_ESP_I3xy2z_F3x+ABZ*I_ESP_H3xyz_F3x;
    abcd[iGrid*315+47] = I_ESP_I3x3z_F3x+ABZ*I_ESP_H3x2z_F3x;
    abcd[iGrid*315+48] = I_ESP_I2x3yz_F3x+ABZ*I_ESP_H2x3y_F3x;
    abcd[iGrid*315+49] = I_ESP_I2x2y2z_F3x+ABZ*I_ESP_H2x2yz_F3x;
    abcd[iGrid*315+50] = I_ESP_I2xy3z_F3x+ABZ*I_ESP_H2xy2z_F3x;
    abcd[iGrid*315+51] = I_ESP_I2x4z_F3x+ABZ*I_ESP_H2x3z_F3x;
    abcd[iGrid*315+52] = I_ESP_Ix4yz_F3x+ABZ*I_ESP_Hx4y_F3x;
    abcd[iGrid*315+53] = I_ESP_Ix3y2z_F3x+ABZ*I_ESP_Hx3yz_F3x;
    abcd[iGrid*315+54] = I_ESP_Ix2y3z_F3x+ABZ*I_ESP_Hx2y2z_F3x;
    abcd[iGrid*315+55] = I_ESP_Ixy4z_F3x+ABZ*I_ESP_Hxy3z_F3x;
    abcd[iGrid*315+56] = I_ESP_Ix5z_F3x+ABZ*I_ESP_Hx4z_F3x;
    abcd[iGrid*315+57] = I_ESP_I5yz_F3x+ABZ*I_ESP_H5y_F3x;
    abcd[iGrid*315+58] = I_ESP_I4y2z_F3x+ABZ*I_ESP_H4yz_F3x;
    abcd[iGrid*315+59] = I_ESP_I3y3z_F3x+ABZ*I_ESP_H3y2z_F3x;
    abcd[iGrid*315+60] = I_ESP_I2y4z_F3x+ABZ*I_ESP_H2y3z_F3x;
    abcd[iGrid*315+61] = I_ESP_Iy5z_F3x+ABZ*I_ESP_Hy4z_F3x;
    abcd[iGrid*315+62] = I_ESP_I6z_F3x+ABZ*I_ESP_H5z_F3x;
    abcd[iGrid*315+63] = I_ESP_I5xy_F2xy+ABY*I_ESP_H5x_F2xy;
    abcd[iGrid*315+64] = I_ESP_I4x2y_F2xy+ABY*I_ESP_H4xy_F2xy;
    abcd[iGrid*315+65] = I_ESP_I4xyz_F2xy+ABY*I_ESP_H4xz_F2xy;
    abcd[iGrid*315+66] = I_ESP_I3x3y_F2xy+ABY*I_ESP_H3x2y_F2xy;
    abcd[iGrid*315+67] = I_ESP_I3x2yz_F2xy+ABY*I_ESP_H3xyz_F2xy;
    abcd[iGrid*315+68] = I_ESP_I3xy2z_F2xy+ABY*I_ESP_H3x2z_F2xy;
    abcd[iGrid*315+69] = I_ESP_I2x4y_F2xy+ABY*I_ESP_H2x3y_F2xy;
    abcd[iGrid*315+70] = I_ESP_I2x3yz_F2xy+ABY*I_ESP_H2x2yz_F2xy;
    abcd[iGrid*315+71] = I_ESP_I2x2y2z_F2xy+ABY*I_ESP_H2xy2z_F2xy;
    abcd[iGrid*315+72] = I_ESP_I2xy3z_F2xy+ABY*I_ESP_H2x3z_F2xy;
    abcd[iGrid*315+73] = I_ESP_Ix5y_F2xy+ABY*I_ESP_Hx4y_F2xy;
    abcd[iGrid*315+74] = I_ESP_Ix4yz_F2xy+ABY*I_ESP_Hx3yz_F2xy;
    abcd[iGrid*315+75] = I_ESP_Ix3y2z_F2xy+ABY*I_ESP_Hx2y2z_F2xy;
    abcd[iGrid*315+76] = I_ESP_Ix2y3z_F2xy+ABY*I_ESP_Hxy3z_F2xy;
    abcd[iGrid*315+77] = I_ESP_Ixy4z_F2xy+ABY*I_ESP_Hx4z_F2xy;
    abcd[iGrid*315+78] = I_ESP_I6y_F2xy+ABY*I_ESP_H5y_F2xy;
    abcd[iGrid*315+79] = I_ESP_I5yz_F2xy+ABY*I_ESP_H4yz_F2xy;
    abcd[iGrid*315+80] = I_ESP_I4y2z_F2xy+ABY*I_ESP_H3y2z_F2xy;
    abcd[iGrid*315+81] = I_ESP_I3y3z_F2xy+ABY*I_ESP_H2y3z_F2xy;
    abcd[iGrid*315+82] = I_ESP_I2y4z_F2xy+ABY*I_ESP_Hy4z_F2xy;
    abcd[iGrid*315+83] = I_ESP_Iy5z_F2xy+ABY*I_ESP_H5z_F2xy;
    abcd[iGrid*315+84] = I_ESP_I5xz_F2xy+ABZ*I_ESP_H5x_F2xy;
    abcd[iGrid*315+85] = I_ESP_I4xyz_F2xy+ABZ*I_ESP_H4xy_F2xy;
    abcd[iGrid*315+86] = I_ESP_I4x2z_F2xy+ABZ*I_ESP_H4xz_F2xy;
    abcd[iGrid*315+87] = I_ESP_I3x2yz_F2xy+ABZ*I_ESP_H3x2y_F2xy;
    abcd[iGrid*315+88] = I_ESP_I3xy2z_F2xy+ABZ*I_ESP_H3xyz_F2xy;
    abcd[iGrid*315+89] = I_ESP_I3x3z_F2xy+ABZ*I_ESP_H3x2z_F2xy;
    abcd[iGrid*315+90] = I_ESP_I2x3yz_F2xy+ABZ*I_ESP_H2x3y_F2xy;
    abcd[iGrid*315+91] = I_ESP_I2x2y2z_F2xy+ABZ*I_ESP_H2x2yz_F2xy;
    abcd[iGrid*315+92] = I_ESP_I2xy3z_F2xy+ABZ*I_ESP_H2xy2z_F2xy;
    abcd[iGrid*315+93] = I_ESP_I2x4z_F2xy+ABZ*I_ESP_H2x3z_F2xy;
    abcd[iGrid*315+94] = I_ESP_Ix4yz_F2xy+ABZ*I_ESP_Hx4y_F2xy;
    abcd[iGrid*315+95] = I_ESP_Ix3y2z_F2xy+ABZ*I_ESP_Hx3yz_F2xy;
    abcd[iGrid*315+96] = I_ESP_Ix2y3z_F2xy+ABZ*I_ESP_Hx2y2z_F2xy;
    abcd[iGrid*315+97] = I_ESP_Ixy4z_F2xy+ABZ*I_ESP_Hxy3z_F2xy;
    abcd[iGrid*315+98] = I_ESP_Ix5z_F2xy+ABZ*I_ESP_Hx4z_F2xy;
    abcd[iGrid*315+99] = I_ESP_I5yz_F2xy+ABZ*I_ESP_H5y_F2xy;
    abcd[iGrid*315+100] = I_ESP_I4y2z_F2xy+ABZ*I_ESP_H4yz_F2xy;
    abcd[iGrid*315+101] = I_ESP_I3y3z_F2xy+ABZ*I_ESP_H3y2z_F2xy;
    abcd[iGrid*315+102] = I_ESP_I2y4z_F2xy+ABZ*I_ESP_H2y3z_F2xy;
    abcd[iGrid*315+103] = I_ESP_Iy5z_F2xy+ABZ*I_ESP_Hy4z_F2xy;
    abcd[iGrid*315+104] = I_ESP_I6z_F2xy+ABZ*I_ESP_H5z_F2xy;
    abcd[iGrid*315+105] = I_ESP_I5xz_F2xz+ABZ*I_ESP_H5x_F2xz;
    abcd[iGrid*315+106] = I_ESP_I4xyz_F2xz+ABZ*I_ESP_H4xy_F2xz;
    abcd[iGrid*315+107] = I_ESP_I4x2z_F2xz+ABZ*I_ESP_H4xz_F2xz;
    abcd[iGrid*315+108] = I_ESP_I3x2yz_F2xz+ABZ*I_ESP_H3x2y_F2xz;
    abcd[iGrid*315+109] = I_ESP_I3xy2z_F2xz+ABZ*I_ESP_H3xyz_F2xz;
    abcd[iGrid*315+110] = I_ESP_I3x3z_F2xz+ABZ*I_ESP_H3x2z_F2xz;
    abcd[iGrid*315+111] = I_ESP_I2x3yz_F2xz+ABZ*I_ESP_H2x3y_F2xz;
    abcd[iGrid*315+112] = I_ESP_I2x2y2z_F2xz+ABZ*I_ESP_H2x2yz_F2xz;
    abcd[iGrid*315+113] = I_ESP_I2xy3z_F2xz+ABZ*I_ESP_H2xy2z_F2xz;
    abcd[iGrid*315+114] = I_ESP_I2x4z_F2xz+ABZ*I_ESP_H2x3z_F2xz;
    abcd[iGrid*315+115] = I_ESP_Ix4yz_F2xz+ABZ*I_ESP_Hx4y_F2xz;
    abcd[iGrid*315+116] = I_ESP_Ix3y2z_F2xz+ABZ*I_ESP_Hx3yz_F2xz;
    abcd[iGrid*315+117] = I_ESP_Ix2y3z_F2xz+ABZ*I_ESP_Hx2y2z_F2xz;
    abcd[iGrid*315+118] = I_ESP_Ixy4z_F2xz+ABZ*I_ESP_Hxy3z_F2xz;
    abcd[iGrid*315+119] = I_ESP_Ix5z_F2xz+ABZ*I_ESP_Hx4z_F2xz;
    abcd[iGrid*315+120] = I_ESP_I5yz_F2xz+ABZ*I_ESP_H5y_F2xz;
    abcd[iGrid*315+121] = I_ESP_I4y2z_F2xz+ABZ*I_ESP_H4yz_F2xz;
    abcd[iGrid*315+122] = I_ESP_I3y3z_F2xz+ABZ*I_ESP_H3y2z_F2xz;
    abcd[iGrid*315+123] = I_ESP_I2y4z_F2xz+ABZ*I_ESP_H2y3z_F2xz;
    abcd[iGrid*315+124] = I_ESP_Iy5z_F2xz+ABZ*I_ESP_Hy4z_F2xz;
    abcd[iGrid*315+125] = I_ESP_I6z_F2xz+ABZ*I_ESP_H5z_F2xz;
    abcd[iGrid*315+126] = I_ESP_I6x_F3y+ABX*I_ESP_H5x_F3y;
    abcd[iGrid*315+127] = I_ESP_I5xy_F3y+ABX*I_ESP_H4xy_F3y;
    abcd[iGrid*315+128] = I_ESP_I5xz_F3y+ABX*I_ESP_H4xz_F3y;
    abcd[iGrid*315+129] = I_ESP_I4x2y_F3y+ABX*I_ESP_H3x2y_F3y;
    abcd[iGrid*315+130] = I_ESP_I4xyz_F3y+ABX*I_ESP_H3xyz_F3y;
    abcd[iGrid*315+131] = I_ESP_I4x2z_F3y+ABX*I_ESP_H3x2z_F3y;
    abcd[iGrid*315+132] = I_ESP_I3x3y_F3y+ABX*I_ESP_H2x3y_F3y;
    abcd[iGrid*315+133] = I_ESP_I3x2yz_F3y+ABX*I_ESP_H2x2yz_F3y;
    abcd[iGrid*315+134] = I_ESP_I3xy2z_F3y+ABX*I_ESP_H2xy2z_F3y;
    abcd[iGrid*315+135] = I_ESP_I3x3z_F3y+ABX*I_ESP_H2x3z_F3y;
    abcd[iGrid*315+136] = I_ESP_I2x4y_F3y+ABX*I_ESP_Hx4y_F3y;
    abcd[iGrid*315+137] = I_ESP_I2x3yz_F3y+ABX*I_ESP_Hx3yz_F3y;
    abcd[iGrid*315+138] = I_ESP_I2x2y2z_F3y+ABX*I_ESP_Hx2y2z_F3y;
    abcd[iGrid*315+139] = I_ESP_I2xy3z_F3y+ABX*I_ESP_Hxy3z_F3y;
    abcd[iGrid*315+140] = I_ESP_I2x4z_F3y+ABX*I_ESP_Hx4z_F3y;
    abcd[iGrid*315+141] = I_ESP_Ix5y_F3y+ABX*I_ESP_H5y_F3y;
    abcd[iGrid*315+142] = I_ESP_Ix4yz_F3y+ABX*I_ESP_H4yz_F3y;
    abcd[iGrid*315+143] = I_ESP_Ix3y2z_F3y+ABX*I_ESP_H3y2z_F3y;
    abcd[iGrid*315+144] = I_ESP_Ix2y3z_F3y+ABX*I_ESP_H2y3z_F3y;
    abcd[iGrid*315+145] = I_ESP_Ixy4z_F3y+ABX*I_ESP_Hy4z_F3y;
    abcd[iGrid*315+146] = I_ESP_Ix5z_F3y+ABX*I_ESP_H5z_F3y;
    abcd[iGrid*315+147] = I_ESP_I5xz_Fx2y+ABZ*I_ESP_H5x_Fx2y;
    abcd[iGrid*315+148] = I_ESP_I4xyz_Fx2y+ABZ*I_ESP_H4xy_Fx2y;
    abcd[iGrid*315+149] = I_ESP_I4x2z_Fx2y+ABZ*I_ESP_H4xz_Fx2y;
    abcd[iGrid*315+150] = I_ESP_I3x2yz_Fx2y+ABZ*I_ESP_H3x2y_Fx2y;
    abcd[iGrid*315+151] = I_ESP_I3xy2z_Fx2y+ABZ*I_ESP_H3xyz_Fx2y;
    abcd[iGrid*315+152] = I_ESP_I3x3z_Fx2y+ABZ*I_ESP_H3x2z_Fx2y;
    abcd[iGrid*315+153] = I_ESP_I2x3yz_Fx2y+ABZ*I_ESP_H2x3y_Fx2y;
    abcd[iGrid*315+154] = I_ESP_I2x2y2z_Fx2y+ABZ*I_ESP_H2x2yz_Fx2y;
    abcd[iGrid*315+155] = I_ESP_I2xy3z_Fx2y+ABZ*I_ESP_H2xy2z_Fx2y;
    abcd[iGrid*315+156] = I_ESP_I2x4z_Fx2y+ABZ*I_ESP_H2x3z_Fx2y;
    abcd[iGrid*315+157] = I_ESP_Ix4yz_Fx2y+ABZ*I_ESP_Hx4y_Fx2y;
    abcd[iGrid*315+158] = I_ESP_Ix3y2z_Fx2y+ABZ*I_ESP_Hx3yz_Fx2y;
    abcd[iGrid*315+159] = I_ESP_Ix2y3z_Fx2y+ABZ*I_ESP_Hx2y2z_Fx2y;
    abcd[iGrid*315+160] = I_ESP_Ixy4z_Fx2y+ABZ*I_ESP_Hxy3z_Fx2y;
    abcd[iGrid*315+161] = I_ESP_Ix5z_Fx2y+ABZ*I_ESP_Hx4z_Fx2y;
    abcd[iGrid*315+162] = I_ESP_I5yz_Fx2y+ABZ*I_ESP_H5y_Fx2y;
    abcd[iGrid*315+163] = I_ESP_I4y2z_Fx2y+ABZ*I_ESP_H4yz_Fx2y;
    abcd[iGrid*315+164] = I_ESP_I3y3z_Fx2y+ABZ*I_ESP_H3y2z_Fx2y;
    abcd[iGrid*315+165] = I_ESP_I2y4z_Fx2y+ABZ*I_ESP_H2y3z_Fx2y;
    abcd[iGrid*315+166] = I_ESP_Iy5z_Fx2y+ABZ*I_ESP_Hy4z_Fx2y;
    abcd[iGrid*315+167] = I_ESP_I6z_Fx2y+ABZ*I_ESP_H5z_Fx2y;
    abcd[iGrid*315+168] = I_ESP_I5xy_Fx2z+ABY*I_ESP_H5x_Fx2z;
    abcd[iGrid*315+169] = I_ESP_I4x2y_Fx2z+ABY*I_ESP_H4xy_Fx2z;
    abcd[iGrid*315+170] = I_ESP_I4xyz_Fx2z+ABY*I_ESP_H4xz_Fx2z;
    abcd[iGrid*315+171] = I_ESP_I3x3y_Fx2z+ABY*I_ESP_H3x2y_Fx2z;
    abcd[iGrid*315+172] = I_ESP_I3x2yz_Fx2z+ABY*I_ESP_H3xyz_Fx2z;
    abcd[iGrid*315+173] = I_ESP_I3xy2z_Fx2z+ABY*I_ESP_H3x2z_Fx2z;
    abcd[iGrid*315+174] = I_ESP_I2x4y_Fx2z+ABY*I_ESP_H2x3y_Fx2z;
    abcd[iGrid*315+175] = I_ESP_I2x3yz_Fx2z+ABY*I_ESP_H2x2yz_Fx2z;
    abcd[iGrid*315+176] = I_ESP_I2x2y2z_Fx2z+ABY*I_ESP_H2xy2z_Fx2z;
    abcd[iGrid*315+177] = I_ESP_I2xy3z_Fx2z+ABY*I_ESP_H2x3z_Fx2z;
    abcd[iGrid*315+178] = I_ESP_Ix5y_Fx2z+ABY*I_ESP_Hx4y_Fx2z;
    abcd[iGrid*315+179] = I_ESP_Ix4yz_Fx2z+ABY*I_ESP_Hx3yz_Fx2z;
    abcd[iGrid*315+180] = I_ESP_Ix3y2z_Fx2z+ABY*I_ESP_Hx2y2z_Fx2z;
    abcd[iGrid*315+181] = I_ESP_Ix2y3z_Fx2z+ABY*I_ESP_Hxy3z_Fx2z;
    abcd[iGrid*315+182] = I_ESP_Ixy4z_Fx2z+ABY*I_ESP_Hx4z_Fx2z;
    abcd[iGrid*315+183] = I_ESP_I6y_Fx2z+ABY*I_ESP_H5y_Fx2z;
    abcd[iGrid*315+184] = I_ESP_I5yz_Fx2z+ABY*I_ESP_H4yz_Fx2z;
    abcd[iGrid*315+185] = I_ESP_I4y2z_Fx2z+ABY*I_ESP_H3y2z_Fx2z;
    abcd[iGrid*315+186] = I_ESP_I3y3z_Fx2z+ABY*I_ESP_H2y3z_Fx2z;
    abcd[iGrid*315+187] = I_ESP_I2y4z_Fx2z+ABY*I_ESP_Hy4z_Fx2z;
    abcd[iGrid*315+188] = I_ESP_Iy5z_Fx2z+ABY*I_ESP_H5z_Fx2z;
    abcd[iGrid*315+189] = I_ESP_I6x_F3z+ABX*I_ESP_H5x_F3z;
    abcd[iGrid*315+190] = I_ESP_I5xy_F3z+ABX*I_ESP_H4xy_F3z;
    abcd[iGrid*315+191] = I_ESP_I5xz_F3z+ABX*I_ESP_H4xz_F3z;
    abcd[iGrid*315+192] = I_ESP_I4x2y_F3z+ABX*I_ESP_H3x2y_F3z;
    abcd[iGrid*315+193] = I_ESP_I4xyz_F3z+ABX*I_ESP_H3xyz_F3z;
    abcd[iGrid*315+194] = I_ESP_I4x2z_F3z+ABX*I_ESP_H3x2z_F3z;
    abcd[iGrid*315+195] = I_ESP_I3x3y_F3z+ABX*I_ESP_H2x3y_F3z;
    abcd[iGrid*315+196] = I_ESP_I3x2yz_F3z+ABX*I_ESP_H2x2yz_F3z;
    abcd[iGrid*315+197] = I_ESP_I3xy2z_F3z+ABX*I_ESP_H2xy2z_F3z;
    abcd[iGrid*315+198] = I_ESP_I3x3z_F3z+ABX*I_ESP_H2x3z_F3z;
    abcd[iGrid*315+199] = I_ESP_I2x4y_F3z+ABX*I_ESP_Hx4y_F3z;
    abcd[iGrid*315+200] = I_ESP_I2x3yz_F3z+ABX*I_ESP_Hx3yz_F3z;
    abcd[iGrid*315+201] = I_ESP_I2x2y2z_F3z+ABX*I_ESP_Hx2y2z_F3z;
    abcd[iGrid*315+202] = I_ESP_I2xy3z_F3z+ABX*I_ESP_Hxy3z_F3z;
    abcd[iGrid*315+203] = I_ESP_I2x4z_F3z+ABX*I_ESP_Hx4z_F3z;
    abcd[iGrid*315+204] = I_ESP_Ix5y_F3z+ABX*I_ESP_H5y_F3z;
    abcd[iGrid*315+205] = I_ESP_Ix4yz_F3z+ABX*I_ESP_H4yz_F3z;
    abcd[iGrid*315+206] = I_ESP_Ix3y2z_F3z+ABX*I_ESP_H3y2z_F3z;
    abcd[iGrid*315+207] = I_ESP_Ix2y3z_F3z+ABX*I_ESP_H2y3z_F3z;
    abcd[iGrid*315+208] = I_ESP_Ixy4z_F3z+ABX*I_ESP_Hy4z_F3z;
    abcd[iGrid*315+209] = I_ESP_Ix5z_F3z+ABX*I_ESP_H5z_F3z;
    abcd[iGrid*315+210] = I_ESP_I5xy_F3y+ABY*I_ESP_H5x_F3y;
    abcd[iGrid*315+211] = I_ESP_I4x2y_F3y+ABY*I_ESP_H4xy_F3y;
    abcd[iGrid*315+212] = I_ESP_I4xyz_F3y+ABY*I_ESP_H4xz_F3y;
    abcd[iGrid*315+213] = I_ESP_I3x3y_F3y+ABY*I_ESP_H3x2y_F3y;
    abcd[iGrid*315+214] = I_ESP_I3x2yz_F3y+ABY*I_ESP_H3xyz_F3y;
    abcd[iGrid*315+215] = I_ESP_I3xy2z_F3y+ABY*I_ESP_H3x2z_F3y;
    abcd[iGrid*315+216] = I_ESP_I2x4y_F3y+ABY*I_ESP_H2x3y_F3y;
    abcd[iGrid*315+217] = I_ESP_I2x3yz_F3y+ABY*I_ESP_H2x2yz_F3y;
    abcd[iGrid*315+218] = I_ESP_I2x2y2z_F3y+ABY*I_ESP_H2xy2z_F3y;
    abcd[iGrid*315+219] = I_ESP_I2xy3z_F3y+ABY*I_ESP_H2x3z_F3y;
    abcd[iGrid*315+220] = I_ESP_Ix5y_F3y+ABY*I_ESP_Hx4y_F3y;
    abcd[iGrid*315+221] = I_ESP_Ix4yz_F3y+ABY*I_ESP_Hx3yz_F3y;
    abcd[iGrid*315+222] = I_ESP_Ix3y2z_F3y+ABY*I_ESP_Hx2y2z_F3y;
    abcd[iGrid*315+223] = I_ESP_Ix2y3z_F3y+ABY*I_ESP_Hxy3z_F3y;
    abcd[iGrid*315+224] = I_ESP_Ixy4z_F3y+ABY*I_ESP_Hx4z_F3y;
    abcd[iGrid*315+225] = I_ESP_I6y_F3y+ABY*I_ESP_H5y_F3y;
    abcd[iGrid*315+226] = I_ESP_I5yz_F3y+ABY*I_ESP_H4yz_F3y;
    abcd[iGrid*315+227] = I_ESP_I4y2z_F3y+ABY*I_ESP_H3y2z_F3y;
    abcd[iGrid*315+228] = I_ESP_I3y3z_F3y+ABY*I_ESP_H2y3z_F3y;
    abcd[iGrid*315+229] = I_ESP_I2y4z_F3y+ABY*I_ESP_Hy4z_F3y;
    abcd[iGrid*315+230] = I_ESP_Iy5z_F3y+ABY*I_ESP_H5z_F3y;
    abcd[iGrid*315+231] = I_ESP_I5xz_F3y+ABZ*I_ESP_H5x_F3y;
    abcd[iGrid*315+232] = I_ESP_I4xyz_F3y+ABZ*I_ESP_H4xy_F3y;
    abcd[iGrid*315+233] = I_ESP_I4x2z_F3y+ABZ*I_ESP_H4xz_F3y;
    abcd[iGrid*315+234] = I_ESP_I3x2yz_F3y+ABZ*I_ESP_H3x2y_F3y;
    abcd[iGrid*315+235] = I_ESP_I3xy2z_F3y+ABZ*I_ESP_H3xyz_F3y;
    abcd[iGrid*315+236] = I_ESP_I3x3z_F3y+ABZ*I_ESP_H3x2z_F3y;
    abcd[iGrid*315+237] = I_ESP_I2x3yz_F3y+ABZ*I_ESP_H2x3y_F3y;
    abcd[iGrid*315+238] = I_ESP_I2x2y2z_F3y+ABZ*I_ESP_H2x2yz_F3y;
    abcd[iGrid*315+239] = I_ESP_I2xy3z_F3y+ABZ*I_ESP_H2xy2z_F3y;
    abcd[iGrid*315+240] = I_ESP_I2x4z_F3y+ABZ*I_ESP_H2x3z_F3y;
    abcd[iGrid*315+241] = I_ESP_Ix4yz_F3y+ABZ*I_ESP_Hx4y_F3y;
    abcd[iGrid*315+242] = I_ESP_Ix3y2z_F3y+ABZ*I_ESP_Hx3yz_F3y;
    abcd[iGrid*315+243] = I_ESP_Ix2y3z_F3y+ABZ*I_ESP_Hx2y2z_F3y;
    abcd[iGrid*315+244] = I_ESP_Ixy4z_F3y+ABZ*I_ESP_Hxy3z_F3y;
    abcd[iGrid*315+245] = I_ESP_Ix5z_F3y+ABZ*I_ESP_Hx4z_F3y;
    abcd[iGrid*315+246] = I_ESP_I5yz_F3y+ABZ*I_ESP_H5y_F3y;
    abcd[iGrid*315+247] = I_ESP_I4y2z_F3y+ABZ*I_ESP_H4yz_F3y;
    abcd[iGrid*315+248] = I_ESP_I3y3z_F3y+ABZ*I_ESP_H3y2z_F3y;
    abcd[iGrid*315+249] = I_ESP_I2y4z_F3y+ABZ*I_ESP_H2y3z_F3y;
    abcd[iGrid*315+250] = I_ESP_Iy5z_F3y+ABZ*I_ESP_Hy4z_F3y;
    abcd[iGrid*315+251] = I_ESP_I6z_F3y+ABZ*I_ESP_H5z_F3y;
    abcd[iGrid*315+252] = I_ESP_I5xz_F2yz+ABZ*I_ESP_H5x_F2yz;
    abcd[iGrid*315+253] = I_ESP_I4xyz_F2yz+ABZ*I_ESP_H4xy_F2yz;
    abcd[iGrid*315+254] = I_ESP_I4x2z_F2yz+ABZ*I_ESP_H4xz_F2yz;
    abcd[iGrid*315+255] = I_ESP_I3x2yz_F2yz+ABZ*I_ESP_H3x2y_F2yz;
    abcd[iGrid*315+256] = I_ESP_I3xy2z_F2yz+ABZ*I_ESP_H3xyz_F2yz;
    abcd[iGrid*315+257] = I_ESP_I3x3z_F2yz+ABZ*I_ESP_H3x2z_F2yz;
    abcd[iGrid*315+258] = I_ESP_I2x3yz_F2yz+ABZ*I_ESP_H2x3y_F2yz;
    abcd[iGrid*315+259] = I_ESP_I2x2y2z_F2yz+ABZ*I_ESP_H2x2yz_F2yz;
    abcd[iGrid*315+260] = I_ESP_I2xy3z_F2yz+ABZ*I_ESP_H2xy2z_F2yz;
    abcd[iGrid*315+261] = I_ESP_I2x4z_F2yz+ABZ*I_ESP_H2x3z_F2yz;
    abcd[iGrid*315+262] = I_ESP_Ix4yz_F2yz+ABZ*I_ESP_Hx4y_F2yz;
    abcd[iGrid*315+263] = I_ESP_Ix3y2z_F2yz+ABZ*I_ESP_Hx3yz_F2yz;
    abcd[iGrid*315+264] = I_ESP_Ix2y3z_F2yz+ABZ*I_ESP_Hx2y2z_F2yz;
    abcd[iGrid*315+265] = I_ESP_Ixy4z_F2yz+ABZ*I_ESP_Hxy3z_F2yz;
    abcd[iGrid*315+266] = I_ESP_Ix5z_F2yz+ABZ*I_ESP_Hx4z_F2yz;
    abcd[iGrid*315+267] = I_ESP_I5yz_F2yz+ABZ*I_ESP_H5y_F2yz;
    abcd[iGrid*315+268] = I_ESP_I4y2z_F2yz+ABZ*I_ESP_H4yz_F2yz;
    abcd[iGrid*315+269] = I_ESP_I3y3z_F2yz+ABZ*I_ESP_H3y2z_F2yz;
    abcd[iGrid*315+270] = I_ESP_I2y4z_F2yz+ABZ*I_ESP_H2y3z_F2yz;
    abcd[iGrid*315+271] = I_ESP_Iy5z_F2yz+ABZ*I_ESP_Hy4z_F2yz;
    abcd[iGrid*315+272] = I_ESP_I6z_F2yz+ABZ*I_ESP_H5z_F2yz;
    abcd[iGrid*315+273] = I_ESP_I5xy_F3z+ABY*I_ESP_H5x_F3z;
    abcd[iGrid*315+274] = I_ESP_I4x2y_F3z+ABY*I_ESP_H4xy_F3z;
    abcd[iGrid*315+275] = I_ESP_I4xyz_F3z+ABY*I_ESP_H4xz_F3z;
    abcd[iGrid*315+276] = I_ESP_I3x3y_F3z+ABY*I_ESP_H3x2y_F3z;
    abcd[iGrid*315+277] = I_ESP_I3x2yz_F3z+ABY*I_ESP_H3xyz_F3z;
    abcd[iGrid*315+278] = I_ESP_I3xy2z_F3z+ABY*I_ESP_H3x2z_F3z;
    abcd[iGrid*315+279] = I_ESP_I2x4y_F3z+ABY*I_ESP_H2x3y_F3z;
    abcd[iGrid*315+280] = I_ESP_I2x3yz_F3z+ABY*I_ESP_H2x2yz_F3z;
    abcd[iGrid*315+281] = I_ESP_I2x2y2z_F3z+ABY*I_ESP_H2xy2z_F3z;
    abcd[iGrid*315+282] = I_ESP_I2xy3z_F3z+ABY*I_ESP_H2x3z_F3z;
    abcd[iGrid*315+283] = I_ESP_Ix5y_F3z+ABY*I_ESP_Hx4y_F3z;
    abcd[iGrid*315+284] = I_ESP_Ix4yz_F3z+ABY*I_ESP_Hx3yz_F3z;
    abcd[iGrid*315+285] = I_ESP_Ix3y2z_F3z+ABY*I_ESP_Hx2y2z_F3z;
    abcd[iGrid*315+286] = I_ESP_Ix2y3z_F3z+ABY*I_ESP_Hxy3z_F3z;
    abcd[iGrid*315+287] = I_ESP_Ixy4z_F3z+ABY*I_ESP_Hx4z_F3z;
    abcd[iGrid*315+288] = I_ESP_I6y_F3z+ABY*I_ESP_H5y_F3z;
    abcd[iGrid*315+289] = I_ESP_I5yz_F3z+ABY*I_ESP_H4yz_F3z;
    abcd[iGrid*315+290] = I_ESP_I4y2z_F3z+ABY*I_ESP_H3y2z_F3z;
    abcd[iGrid*315+291] = I_ESP_I3y3z_F3z+ABY*I_ESP_H2y3z_F3z;
    abcd[iGrid*315+292] = I_ESP_I2y4z_F3z+ABY*I_ESP_Hy4z_F3z;
    abcd[iGrid*315+293] = I_ESP_Iy5z_F3z+ABY*I_ESP_H5z_F3z;
    abcd[iGrid*315+294] = I_ESP_I5xz_F3z+ABZ*I_ESP_H5x_F3z;
    abcd[iGrid*315+295] = I_ESP_I4xyz_F3z+ABZ*I_ESP_H4xy_F3z;
    abcd[iGrid*315+296] = I_ESP_I4x2z_F3z+ABZ*I_ESP_H4xz_F3z;
    abcd[iGrid*315+297] = I_ESP_I3x2yz_F3z+ABZ*I_ESP_H3x2y_F3z;
    abcd[iGrid*315+298] = I_ESP_I3xy2z_F3z+ABZ*I_ESP_H3xyz_F3z;
    abcd[iGrid*315+299] = I_ESP_I3x3z_F3z+ABZ*I_ESP_H3x2z_F3z;
    abcd[iGrid*315+300] = I_ESP_I2x3yz_F3z+ABZ*I_ESP_H2x3y_F3z;
    abcd[iGrid*315+301] = I_ESP_I2x2y2z_F3z+ABZ*I_ESP_H2x2yz_F3z;
    abcd[iGrid*315+302] = I_ESP_I2xy3z_F3z+ABZ*I_ESP_H2xy2z_F3z;
    abcd[iGrid*315+303] = I_ESP_I2x4z_F3z+ABZ*I_ESP_H2x3z_F3z;
    abcd[iGrid*315+304] = I_ESP_Ix4yz_F3z+ABZ*I_ESP_Hx4y_F3z;
    abcd[iGrid*315+305] = I_ESP_Ix3y2z_F3z+ABZ*I_ESP_Hx3yz_F3z;
    abcd[iGrid*315+306] = I_ESP_Ix2y3z_F3z+ABZ*I_ESP_Hx2y2z_F3z;
    abcd[iGrid*315+307] = I_ESP_Ixy4z_F3z+ABZ*I_ESP_Hxy3z_F3z;
    abcd[iGrid*315+308] = I_ESP_Ix5z_F3z+ABZ*I_ESP_Hx4z_F3z;
    abcd[iGrid*315+309] = I_ESP_I5yz_F3z+ABZ*I_ESP_H5y_F3z;
    abcd[iGrid*315+310] = I_ESP_I4y2z_F3z+ABZ*I_ESP_H4yz_F3z;
    abcd[iGrid*315+311] = I_ESP_I3y3z_F3z+ABZ*I_ESP_H3y2z_F3z;
    abcd[iGrid*315+312] = I_ESP_I2y4z_F3z+ABZ*I_ESP_H2y3z_F3z;
    abcd[iGrid*315+313] = I_ESP_Iy5z_F3z+ABZ*I_ESP_Hy4z_F3z;
    abcd[iGrid*315+314] = I_ESP_I6z_F3z+ABZ*I_ESP_H5z_F3z;
  }
}
