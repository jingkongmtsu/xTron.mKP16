#ifndef SPHEREINT
#define SPHEREINT
#include<vector>
#include "libgen.h"

using namespace std;

namespace sphereint
{
	const UInt DIM(3);

	void sphereInt(const LInt& Lcode, const UInt& inp2, const UInt& nGrids, 
	     const Double* icoe, const Double* iexp, const Double* P, const Double* A, 
	     const Double* B, const Double* R, Double* abcd, Double sValue);

	Double gaussian(Int i, Int j, Int k, Double rho, Double norms, vector<Double>& p, 
	       vector<Double>& r);
	Double regressive(Int i, Int j, Int k, Double mu, Double normq, Double norms, 
	       vector< vector<Double> >& invA);
	void getexpanM(Int i, Int j, Int k, Double normq, vector< vector<Double> >& invA, 
	     vector< vector< vector<Double> > >& expanM);
	void getInvM(Double theta0, Double phi0, vector< vector<Double> >& invA);
	Double expanFun(Int i, Int j, Int k, vector<Int> position, 
	       vector< vector<Double> >& invA, Int nVariable, Double normq);
	Double integralF(const Int& i, const Int& j, const Int& k, const Int& l, 
	       const Double& u);
	Double intThetaOdd(const Int& m, const Double& mu);
	//Double intThetaEven(const Int& m, const Double& mu);
	Double approchBessel0(const Double& mu);
	Double approchBessel1(const Double& mu);
	Double approchBessel2(const Double& mu);
	Double expanFun(const Int& n, const Double& mu);
	Double erf(const Double& mu);
	//Double intPhi(const Int& n);
	Int factorial(const Int& n);

	//Calculate the factor Fijk showed in fenlai's note formula 1.69
	Double coefCalculator(Int& i, Int& j, Int& k, Int& l1, Int& l2, Int& m1, Int& m2, 
	       Int& n1, Int& n2, vector<Double> PA, vector<Double> PB);
	Double coefCal(Int k, Int l1, Int l2, Double PA, Double PB);
	Int combination(Int n, Int k);
}

#endif 
