/**
 * \file   functions.h
 * \brief  This file contains the variety mathematical functions 
 */
#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "general.h"
#include <boost/math/special_functions/binomial.hpp>  // special functions in math
#include <boost/math/special_functions/factorials.hpp>
using namespace boost::math;

namespace functions {

	/**
	 * factorial of n!
	 */
	inline Double getFactorial(const Double& n) {
		return factorial<Double>(n);
	};

	/**
	 * BionomialCoe returns C^{l}_{n}
	 */
	inline Double bionomialCoe(const Double& n, const Double& l) {
		return binomial_coefficient<Double>(n, l);
	};

	/**
	 * DoubleFac returns N!!
	 */
	inline Double doubleFac(const Double& n) {
		if (fabs(n)<THRESHOLD_MATH || fabs(n+ONE)<THRESHOLD_MATH) {
			return ONE; //case for n = -1 and n = 0
		}else{
			return double_factorial<Double>(n);
		}
	}; 


}

#endif
