#include<stdlib.h>
#include<iostream>
#include<stdio.h>
#include<cmath>
#include<vector>

using namespace std;

//#define Double double
//#define Int int

#include"sphereint.h"

using std::max;
using std::min;

using namespace sphereint;

// This is a function to calculate factor Fijk which showed in
// fenlai's note in formula 1.69


Double sphereint::coefCalculator(Int& i, Int& j, Int& k, Int& l1, Int& l2, Int& m1, Int& m2, Int& n1, Int& n2, vector<Double> PA, vector<Double> PB){
  Double result = Double(1);

  Double xPA = PA[0];
  Double yPA = PA[1];
  Double zPA = PA[2];
  Double xPB = PB[0];
  Double yPB = PB[1];
  Double zPB = PB[2];
  
  result *= coefCal(i, l1, l2, xPA, xPB);

  result *= coefCal(j, m1, m2, yPA, yPB);

  result *= coefCal(k, n1, n2, zPA, zPB);
  
  return result;

}
// The formula from fenlai's note

Double sphereint::coefCal(Int k, Int l1, Int l2, Double PA, Double PB){
  Double result = 0;
  if(l1==0)
    {
      if(l2==0)
	{
	  result = 1;
	  return result;
	}
      else
	{
	  result = pow(PB,l2-k)*combination(l2,k);
	  return result;
	}
    }
  else
    {
      if(l2==0)
	{
	  result = pow(PA,l1-k)*combination(l1,k);
	  return result;
	}
      else
	{
	  Int upper = min(k,l1);
	  Int lowwer = max(0,k-l2);

	  for(Int i=lowwer; i<upper+1; i++)
	    {
	      result += pow(PA,l1-i)*pow(PB,l2-k+i)*combination(l1,i)*combination(l2,k-i);
	    }

	  return result;
	}
    }
}

// This is combination function  C(n,k) = n!/(k! * (n-k)!)
Int sphereint::combination(Int n, Int k){
  Int result1 = 1;
  Int result2 = 1;

  if(k > n   || k < 0){
    cout << "error entry" << endl;
    return 1;
  }

  for(Int i=1; i<k+1; i++){
    result1 *= n-i+1;
    result2 *= i;
  }
  
  return result1/result2;

}
