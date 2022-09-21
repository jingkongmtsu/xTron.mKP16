/**
 * \file   functions.h
 * \brief  This file contains the variety mathematical functions 
 * \author Fenglai Liu 
 */
#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "libgen.h"
#include <boost/math/special_functions/binomial.hpp>  // special functions in math
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace functions {

	/**
	 * factorial of n!
	 */
	inline Double getFactorial(Double n) {
		return boost::math::factorial<Double>(n);
	};

	/**
	 * BionomialCoe returns C^{l}_{n}
	 */
	inline Double bionomialCoe(Double n, Double l) {
		return boost::math::binomial_coefficient<Double>(n, l);
	};

	/**
	 * DoubleFac returns N!!
	 */
	inline Double doubleFac(Double n) {
		if (fabs(n)<THRESHOLD_MATH || fabs(n+ONE)<THRESHOLD_MATH) {
			return ONE; //case for n = -1 and n = 0
		}else{
			return boost::math::double_factorial<Double>(n);
		}
	}; 

	/**
	 * GammaFunInt calculate the Gamma integral below:
	 * \int_-infinity^+infinity  x^n*exp(-alpha*x^2) dx
	 *
	 * Be careful about the case that alpha is near zero
	 * here we did not check it
	 */
   inline Double gammaFunInt(UInt n, Double alpha) {
		if (n%2 == 0) {
			if (n == 0) {
				return sqrt(PI/alpha);
			} else {
				Double x  = sqrt(PI/alpha);
				Double ta = TWO*alpha;
				for(UInt i=n/2; i>0; i--) {
					x = x*(TWO*i-1)/ta;
				}
				return x;
			}
		} else {
			return ZERO;
		}
	};

}

#endif
