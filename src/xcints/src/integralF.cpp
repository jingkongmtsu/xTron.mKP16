// This function is used to calculate the integral of the following function
// F(theta,phi) = integral[(sin(theta)^i)*(cos(theta)^j)*(sin(phi)^k)*(cos(phi)^l)*exp(-mu*cos(theta)),{theta,0,Pi},{phi,0,2*Pi}]
// There are 5 input variable:
// Int: i,j,k,l
// Double : mu
// The output is Double

#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<cmath>
#include<vector>
#include<gsl/gsl_sf_bessel.h>
using namespace std;
using std::max;
using std::min;
#define BESSELBOND 500  // (* This one may change a little bit smaller if the number is overlap)

#include"sphereint.h"
#include"shellprop.h"

using namespace sphereint;


Int sphereint::factorial(const Int& n){
  if(n < 0){
    cout << "error, n is less than 0" << endl;
    return 0;
  }
  if(n == 0 || n == 1){
    return 1;
  }else{
    return n*sphereint::factorial(n-1);
  }
}


Double intPhi(const Int& n){
  if(n % 2 == 1){
    return 0;
  }else if(n == 0){
    return 4*PI;
  }else{
    return (n-1)*1.0/n*intPhi(n-2);
  }
}

Double intThetaEven(const Int& m, const Double& mu){
  if(m > 2){
    return ((m-1)*intThetaEven(m-1,mu) + mu*intThetaEven(m-2,mu) - (m-2)*intThetaEven(m-3,mu)) / mu;
  }else if(m == 2){
    Double ans2;
    if(mu < BESSELBOND){
      ans2 = gsl_sf_bessel_In(1,mu)/exp(mu)*PI/mu + gsl_sf_bessel_In(2,mu)/exp(mu)*PI;
    }else{
      ans2 = approchBessel2(mu);
    }
    return ans2;
  }else if(m == 1){
    Double ans1;
    if(mu < BESSELBOND){
      ans1 = -gsl_sf_bessel_In(1,mu)/exp(mu)*PI;
    }else{
      ans1 = approchBessel1(mu);
    }
    return ans1;
  }else{
    Double ans0;
    if(mu < BESSELBOND){
      ans0 = gsl_sf_bessel_In(0,mu)/exp(mu)*PI;
    }else{
      ans0 = approchBessel0(mu);
    }
    return ans0;
  }    
}

Double sphereint::integralF(const Int& i, const Int& j, const Int& k, const Int& l, const Double& mu){
  // Test if i,j,k,l >= 0. If not, return error message.
  if(i < 0 || j < 0 || k < 0 || l < 0){
    cout << "There is at least one negative value in i,j,k,l. Please check again!" << endl;
    return 0;
  }
  // If all i,j,k,l greater than 0, continue do the work
  else{

    if(k % 2 == 1 || l % 2 == 1){
      return 0;  // The value for intPhi = 0, see formula 1.27

    }else{
      Double H = Double(0.0); // The value for intTheta function
      if(i % 2 == 1){
	// This part calculate the H function when i is an odd number
	Int i1 = (i-1)/2;
	vector<Int> fij; // Define a vector to save the coefficient for H
	
	fij.reserve(i1+1);
	fij.push_back(-1);

	// Construct the coef for fij
	for(Int s=0; s<i1; s++){
	  Double tempfij;
	  tempfij = (i1-s)*fij[s]/(s+1);
	  tempfij = (-1)*tempfij;
	  fij.push_back(tempfij);
	}

	// Accumulate the value for H
	for(Int s=0; s<i1+1; s++){
	  H += fij[s]*intThetaOdd(j+2*s,mu);
	}	
      }else{ 
	// This part calculate the H function when i is an even number
	Int i1 = i/2;

	vector<Int> fij; // vector to save the coefficient for H

	fij.reserve(i1+1);
	fij.push_back(1);

	// Construct the coef for fij
	for(Int s=0; s<i1; s++){
	  Double tempfij;
	  tempfij = (i1-s)*fij[s]/(s+1);
	  tempfij = (-1)*tempfij;
	  //cout << "s = " << s << " ,factor = " << tempfij << endl;
	  fij.push_back(tempfij);
	}

	// Accumulate the value for H
	for(Int s=0; s<i1+1; s++){
	  H += fij[s]*intThetaEven(j+2*s,mu);
	}	
      }

      Double G = Double(0.0); // The value for intPhi function

      Int k1 = k/2;
      Int l1 = l/2;

      vector<Int> gkl; // vector to save the coefficient for G

      gkl.reserve(k1+l1+1);

      Double temgkl = Double(1.0)/pow(Double(2),k1+l1+1); 
      // Construct the coef for gkl
      for(Int t=0; t<k1+l1+1;t++){
	Double tempgkl = Double(0);
	Int uperBond = min(k1,t);
	Int lowerBond = max(0,t-l1);
	for(Int s=lowerBond; s<uperBond+1; ++s){
	  tempgkl += pow(-1,s)*factorial(k1)*factorial(l1)/(factorial(s)*factorial(k1-s))/(factorial(t-s)*factorial(l1-t+s));
	}
	gkl.push_back(tempgkl);
      }
      
      // Accumulate the value for G
      for(Int t=0; t<k1+l1+1; t++){
	G += gkl[t]*intPhi(t);
      }

      G *= temgkl;

      return H*G;
    }
  }
}

Double sphereint::intThetaOdd(const Int& m, const Double& mu){
  vector<Double> coef;
  
  coef.reserve(m+1);

  Double temCoef;

  temCoef = factorial(m)*Double(1.0)/pow(mu,m+1);
  coef.push_back(temCoef);

  for(Int i=1; i<m+1; i++){
    temCoef = coef[i-1]/i*mu;
    coef.push_back(temCoef);    
  }

  Double int1 = Double(0);
  Double int2 = Double(0);
  for(Int i=0; i<m+1; i++){
    int1 += coef[i];
    int2 += coef[i]*pow(-1,i);
  }

  int1 *= exp(-Double(2)*mu);
  return int1 - int2;
}
      
Double sphereint::approchBessel0(const Double& mu){
  Int expanNum = 5; // the number of expansion;
  Double ans = Double(0.0);

  Double factor;
  for(Int n=0; n<expanNum; n++){
    Double tempAns = expanFun(n,mu);
    factor = 1;
    for(Int m=1; m<n+1; ++m){
      factor = factor*(Double(2)*m-1);
    }
    factor = factor*Double(1.0)/pow(Double(2),n)/factorial(n);
    ans = ans + tempAns*factor;
  }
  ans = ans*2;
  return ans;
}

Double sphereint::erf(const Double& mu){
  Double ans;
  ans = sqrt(PI)/Double(2) + Double(31)/Double(200)*exp(-Double(2)*mu) - Double(341)/Double(8000)*exp(-Double(2)*Double(2)*mu);
  ans *= Double(2)/sqrt(PI)*sqrt(1-exp(-Double(2)*mu));
  return ans;
}

Double sphereint::expanFun(const Int& n, const Double& mu){
  Double expan = Double(0.0);
  switch(n){
  case 0:
    expan = sqrt(PI/2)*erf(mu)*Double(0.5)/sqrt(mu);
    break;
  case 1:
    expan = (-4*exp(-2*mu)*sqrt(mu) + sqrt(2*PI)*erf(mu)) / (16*pow(mu,1.5));
    break;
  case 2:
    expan = (-4*exp(-2*mu)*sqrt(mu)*(Double(3)+4*mu) + 3*sqrt(2*PI)*erf(mu)) / (64*pow(mu,2.5));
    break;
  case 3:
    expan = (-4*exp(-2*mu)*sqrt(mu)*(Double(15)+4*mu*(5+4*mu)) + 15*sqrt(2*PI)*erf(mu)) / (256*pow(mu,3.5));
    break;
  case 4:
    expan = (-4*exp(-2*mu)*sqrt(mu)*(Double(105)+4*mu*(35+4*mu*(7+4*mu))) + 105*sqrt(2*PI)*erf(mu)) / (1024*pow(mu,4.5));
    break;
  default:
    cout << "Expansion exceed the expansion number." << endl;
    break;
  }
  return expan;
}

Double sphereint::approchBessel1(const Double& mu){
  Int expanNum = 5; // the number of expansion;
  Double ans = Double(0.0);

  Double factor, factor1, factor2;
  for(Int n=0; n<expanNum; n++){
    Double tempAns = expanFun(n,mu);
    if(n == 0){
      factor1 = 0;
    }else{
      factor1 = 1;
      for(Int m=1; m<n; ++m){
	factor1 = factor1*(2*m-1);
      }
      factor1 = 2*factor1*1.0/pow(2,n-1)/factorial(n-1);
    }
    
    factor2 = 1;
    for(Int m=1; m<n+1; ++m){
      factor2 = factor2*(2*m-1);
    }
    factor2 = factor2*1.0/pow(2,n)/factorial(n);
    
    factor = factor1 - factor2;
    ans = ans + tempAns*factor;
  }
  ans = ans*2;
  return ans;
}

Double sphereint::approchBessel2(const Double& mu){
  Int expanNum = 5; // the number of expansion;
  Double ans = 0.0;

  Double factor, factor1, factor2, factor3;
  for(Int n=0; n<expanNum; n++){
    Double tempAns = expanFun(n,mu);

    if(n == 0){
      factor1 = 0;
    }else if(n == 1){
      factor1 = 0;
    }else{
      factor1 = 1;
      for(Int m=1; m<n-1; ++m){
	factor1 = factor1*(2*m-1);
      }
      factor1 = 4*factor1*1.0/pow(2,n-2)/factorial(n-2);
    }

    if(n == 0){
      factor2 = 0;
    }else{
      factor2 = 1;
      for(Int m=1; m<n; ++m){
	factor2 = factor2*(2*m-1);
      }
      factor2 = 4*factor2*1.0/pow(2,n-1)/factorial(n-1);
    }

    factor3 = 1;
    for(Int m=1; m<n+1; ++m){
      factor3 = factor3*(2*m-1);
    }
    factor3 = factor3*1.0/pow(2,n)/factorial(n);

    //cout << "factor1 = " << factor1 << endl;
    //cout << "factor2 = " << factor2 << endl;
    //cout << "factor3 = " << factor3 << endl;

    factor = factor1 - factor2 + factor3;
    
    ans = ans + tempAns*factor;
  }
  
  ans = ans*2;

  return ans;
}

