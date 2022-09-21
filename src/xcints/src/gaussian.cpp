#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<vector>
#include<cmath>

#include "sphereint.h"
//#include "integral.h"

using namespace std;
using std::min;
using std::max;
using namespace sphereint;

Double sphereint::gaussian(Int i, Int j, Int k, Double rho, Double norms, vector<Double>& p, vector<Double>& r){
  Double normq;
  vector<Double> q;  // the vector q = r - p;

  q.reserve(DIM);

  for(Int t=0; t<DIM; t++){
    q.push_back(r[t]-p[t]);
  }

  normq = Double(0.0); // Initialize the norm of vector q

  // Get the mod length for vector s ans q
  for(Int t=0; t<DIM; t++){
    normq += q[t]*q[t];
  }

  normq = sqrt(normq);

  // This part is to get trigonometric functions
  // see formula (1.2)
  Double theta0, phi0; // the theta0, phi0 is the angle vector q
  theta0 = acos(q[2]/normq);
  phi0 = atan2(q[1],q[0]);

  vector< vector<Double> > invA;

  //This function is to get the inverse matrix for A
  //see formula 1.4 and 1.5

  getInvM(theta0, phi0, invA);

  Double result = 0;
  Double gamma;

  gamma = rho;

  Double factor;

  factor = exp(-gamma*(norms*norms + normq*normq) + 2*gamma*norms*normq);

  Double mu;

  mu = 2*gamma*norms*normq;

  result = factor * regressive(i, j, k, mu, normq, norms, invA);
      
  return result;
  
}

// This function is used for get thr inverse matrix to transfer the cartissian
void sphereint::getInvM(Double theta0, Double phi0, vector< vector<Double> >& invA){
  invA.reserve(DIM);
  vector<Double> temp;
  temp.reserve(DIM);

  // input the matrix A
  // see formula 1.4
  temp.push_back(cos(theta0) * cos(phi0));
  temp.push_back(-sin(phi0));
  temp.push_back(sin(theta0) * cos(phi0));
  invA.push_back(temp);
  temp.clear();

  temp.push_back(cos(theta0) * sin(phi0));
  temp.push_back(cos(phi0));
  temp.push_back(sin(theta0) * sin(phi0));
  invA.push_back(temp);
  temp.clear();

  temp.push_back(-sin(theta0));
  temp.push_back(Double(0));
  temp.push_back(cos(theta0));
  invA.push_back(temp);
  temp.clear();
}

// This function is used for accumlate the functions
Double sphereint::regressive(Int i, Int j, Int k, Double mu, Double normq, Double norms, vector< vector<Double> >& invA){
  vector< vector< vector<Double> > > expanM;
  vector< vector<Double> > tempv;
  vector<Double> temp;
  Double result = 0;

  expanM.reserve(i+j+k+1);
  tempv.reserve(i+j+k+1);
  temp.reserve(i+j+k+1);

  //Initial expansion matrix
  for(Int t=0; t<i+j+k+1; t++){
    temp.push_back(Double(0));
  }
  
  for(Int t=0; t<i+j+k+1; t++){
    tempv.push_back(temp);
  }

  for(Int t=0; t<i+j+k+1; t++){
    expanM.push_back(tempv);
  }
  
  getexpanM(i, j, k, normq, invA, expanM);
  
  Int upper = i+j+k+1;  

  for(Int t=0; t<upper; t++){
    for(Int s=0; s<upper-t; ++s){
      for(Int r=0; r<upper-s-t; ++r){
	result += expanM[t][s][r] * pow(norms,t+s+r) * integralF(t+s+1,r,s,t,mu);
      }
    }
  }
  
  return result; 
  
}

// This function try to get the expand function in each position
void sphereint::getexpanM(Int i, Int j, Int k, Double normq, vector< vector<Double> >& invA, vector< vector< vector<Double> > >& expanM){

  Int upper = i+j+k+1;
  Int nVariable = DIM;
  for(Int t=0; t<upper; t++){
    for(Int s=0; s<upper-t; ++s){
      for(Int r=0; r<upper-s-t; ++r){
	vector<Int> position;
	position.reserve(DIM);
	position.push_back(t);
	position.push_back(s);
	position.push_back(r);
	expanM[t][s][r] = expanFun(i,j,k,position,invA,nVariable,normq);
      }
    }
  }
}

// This is the expand function, also the regressive function 
// See formula 1.7
Double sphereint::expanFun(Int i, Int j, Int k, vector<Int> position, vector< vector<Double> >& invA, Int nVariable, Double normq){
  Int upperBondl, upperBondm, lowerBondl, lowerBondm;
  Double factor = 0;

  if(nVariable == 0){
    Double a, b, c;
    a = invA[0][2]*normq;
    b = invA[1][2]*normq;
    c = invA[2][2]*normq;
    
    return pow(a,i) * pow(b,j) * pow(c,k);
  }
  
  Double a, b, c;
  Int t;
  a = invA[0][DIM - nVariable];
  b = invA[1][DIM - nVariable];
  c = invA[2][DIM - nVariable];
  t = position[DIM - nVariable];
  upperBondl = min(i,t) + 1;
  lowerBondl = max(0,t-(j+k));
  for(Int l=lowerBondl; l<upperBondl; ++l){
    upperBondm = min(t-l,j) + 1;
    lowerBondm = max(0,t-l-k);
    for(Int m=lowerBondm; m<upperBondm; ++m){
      Double tempfactor;
      tempfactor = pow(a,l) * pow(b,m) * pow(c,t-l-m);
      tempfactor *= factorial(i)/factorial(l)/factorial(i-l);
      tempfactor *= factorial(j)/factorial(m)/factorial(j-m);
      tempfactor *= factorial(k)/factorial(t-l-m)/factorial(k-t+l+m);
      factor += tempfactor*sphereint::expanFun(i-l, j-m, k+l+m-t, position, invA, nVariable-1, normq);
    }
  }
 
  return factor;
}
