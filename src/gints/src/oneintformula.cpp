/**
 * CPP files corresponding to the one electron integrals calculation 
 * in direct way
 * \author fenglai liu
 */
#include<cmath>
#include<cstdio>
#include<iostream>
#include "functions.h"
#include "oneintformula.h"
using namespace functions;
using namespace oneintformula;

Double oneintformula::Fk(const UInt& k, const UInt& l1, const UInt& l2, 
			const Double& x1, const Double& x2) 
{
	//
	// firstly, it's necessary to note that we can not use UInt inside
	// the calculation. We may produce negative number during the 
	// calculation. Therefore, let's transfer the UInt into Int first
	//
	Int k_  = static_cast<Int>(k);
	Int l1_ = static_cast<Int>(l1);
	Int l2_ = static_cast<Int>(l2);

	//
	// we do not check that whether l1, l2 are >=0
	// this is assumed to be true
	// but you need to be aware of this!
	//
	if (l1_ >0 && l2_ == 0) {
		return pow(x1,l1_-k_)*bionomialCoe(l1_,k_);
	} else if (l1_ ==0 && l2_ > 0) {
		return pow(x2,l2_-k_)*bionomialCoe(l2_,k_);
	} else if (l1_ ==0 && l2_ == 0) {
		return ONE;
	} 

	// now the only possible case is that l1 > 0 and l2 > 0
	//UInt low   = std::max(0,k-l2);
	//UInt upper = std::min(k,l1);
	Int low   = (0<k_-l2_ ? k_-l2_:0);
	Int upper = (k_<l1_   ? k_:l1_);
	Double sum = ZERO;
	for (Int i=low;i<=upper;i++) {
		sum += pow(x1,l1_-i)*pow(x2,l2_-k_+i)*bionomialCoe(l1_,i)*bionomialCoe(l2_,k_-i);
	}
	return sum;
}

Double oneintformula::overlapIntegral(
			const Double& d,   const Double& alpha, 
			const Double& PAx, const Double& PAy, const Double& PAz,  
			const Double& PBx, const Double& PBy, const Double& PBz,  
			const UInt& li, const UInt& mi, const UInt& ni,
			const UInt& lj, const UInt& mj, const UInt& nj, bool sameCenter) 
{
	//
	// for the same center case, we need to treat it here
	// the way below for computing Fk will give all zero value
	// since all PAi and PBi are zero for same center
	//
	if (sameCenter) {
		return d*gammaFunInt(li+lj,alpha)*gammaFunInt(mi+mj,alpha)*gammaFunInt(ni+nj,alpha);
	}

	//
	// Ix
	//
	Double Ix  = ZERO;
	for (UInt i=0; i<=(UInt)(li+lj)/2; i++) {
		Double fk = Fk(2*i,li,lj,PAx,PBx);
		Double g  = gammaFunInt(2*i,alpha);
		Ix += fk*g; 
	}

	//
	// Iy
	//
	Double Iy  = ZERO;
	for (UInt i=0; i<=(UInt)(mi+mj)/2; i++) {
		Double fk= Fk(2*i,mi,mj,PAy,PBy);
		Double g = gammaFunInt(2*i,alpha);
		Iy += fk*g; 
	}

	//
	// Iz
	//
	Double Iz  = ZERO;
	for (UInt i=0; i<=(UInt)(ni+nj)/2; i++) {
		Double fk= Fk(2*i,ni,nj,PAz,PBz);
		Double g = gammaFunInt(2*i,alpha);
		Iz += fk*g; 
	}

	// now it's result
	return d*Ix*Iy*Iz;

}

