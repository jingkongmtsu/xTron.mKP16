#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<vector>


using namespace std;

#include "shellprop.h"
//#include "gaussian.h"
//#include "coefCalculator.h"
#include "sphereint.h"
#include "shell.h"
#include "shellpair.h"

using namespace shellprop;
using namespace sphereint;

#define Double double
#define Int int



void sphereint::sphereInt(const LInt& Lcode, const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd, Double sValue){

  //loop over grid points
  for(UInt iGrid=0; iGrid<nGrids; iGrid++){
    
    // get how many base functions from A and B
    UInt lAmin,lAmax,lBmin,lBmax;
    decodeL(Lcode,lAmin,lAmax,lBmin,lBmax); // Get the smallest and lagest l value

    UInt getBasA = getCartBas(lAmin,lAmax); // Get the number of basic function for A

    UInt getBasB = getCartBas(lBmin,lBmax); // Get the number of basic function for B

    UInt totalBas = getBasA*getBasB; // The total combination for A and B


      for(UInt indexB=0; indexB<getBasB; indexB++){

	UInt lB,mB,nB;
	getlmn(lBmin,lBmax,indexB,lB,mB,nB); // Get the power position
	Int int_lB = lB;
	Int int_mB = mB;
	Int int_nB = nB;

	for(UInt indexA=0; indexA<getBasA; indexA++){
	  UInt lA,mA,nA;
	  getlmn(lAmin,lAmax,indexA,lA,mA,nA); // Get the power position
	  Int int_lA = lA;
	  Int int_mA = mA;
	  Int int_nA = nA;
	  //cout <<"(indexA,lA,mA,nA) = " << indexA<<" "<<lA<<" "<<mA<<" "<<nA<<endl;
	  
	  for(UInt ip2=0; ip2<inp2; ip2++){
	    //Intial number 
	    Double ic2 = icoe[ip2]; // The coefficient 
	    Double onedz = iexp[ip2];
	    Double rho = 1.0E0/onedz; // rho = alpha + beta
	    //Double iexdiff = iexpdiff[ip2];
	    //Double alpha = (rho + iexdiff)/2;
	    //Double beta = (rho + iexdiff)/2;
	    UInt offsetP = 3*ip2;
	    UInt offsetR = 3*iGrid;
	    Double PX = P[offsetP ];
	    Double PY = P[offsetP+1];
	    Double PZ = P[offsetP+2];
	    Double RX = R[offsetR ];
	    Double RY = R[offsetR+1];
	    Double RZ = R[offsetR+2];
	    Double PAX = PX - A[0];
	    Double PAY = PY - A[1];
	    Double PAZ = PZ - A[2];
	    Double PBX = PX - B[0];
	    Double PBY = PY - B[1];
	    Double PBZ = PZ - B[2];
            //Double ABX = A[0] - B[0];
	    //Double ABY = A[1] - B[1];
	    //Double ABZ = A[2] - B[2];
	    //Double AB2 = ABX*ABX + ABY*ABY + ABZ*ABZ;
	    Double prefactor = ic2;
	          
	    vector<Double> p(DIM); // The vector for center of A and B
	    p[0] = PX;
	    p[1] = PY;
	    p[2] = PZ;

	    vector<Double> r(DIM); // The reference point
	    r[0] = RX;
	    r[1] = RY;
	    r[2] = RZ;

	    vector<Double> PA(DIM);
	    vector<Double> PB(DIM);
	    //PA = (Double*)malloc(DIM);
	    //PB = (Double*)malloc(DIM);

	    PA[0] = PAX;
	    PA[1] = PAY;
	    PA[2] = PAZ;
	    PB[0] = PBX;
	    PB[1] = PBY;
	    PB[2] = PBZ;

	    Double result = 0;
	 
	    //Do the accumulation , see formula 1.69 from fenlai's note
	    for(UInt i=0; i<lA+lB+1; i++){
	      for(UInt j=0; j<mA+mB+1; j++){
		for(UInt k=0; k<nA+nB+1; k++){
		  Int int_i = i;
		  Int int_j = j;
		  Int int_k = k;
    
		  result += gaussian(int_i,int_j,int_k,rho,sValue,p,r)*coefCalculator(int_i,int_j,int_k,int_lA,int_lB,int_mA,int_mB,int_nA,int_nB,PA,PB);

		}
	      }
	    }
	
	    // Multiply the prefactor for the result
	    result *= prefactor;

	    abcd[iGrid*totalBas + indexB*getBasA + indexA] += result;

	  }
	}
      }
    
  }
}

