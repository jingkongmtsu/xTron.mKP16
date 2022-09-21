/**
 * \file  oneintformula.h
 * \brief This file contains the math formulas for single electron integral
 * \author Fenglai Liu 
 */
#ifndef ONEINTFORMULA_H
#define ONEINTFORMULA_H
#include "libgen.h"

namespace oneintformula {

	/**
	 * This function is used to calculate the fk in the form like below:
	 *
	 * \f$(X + a)^{l1}(X + b)^{l2} = \sum_{k=0}^{l1+l2}X^{k}fk(a,b,l1,l2)\f$
	 *
	 * This expression is used in the Gaussian Primitive Product Theorem
	 * A, B are atom centers, P is the new center created in the theorem
	 * l are the original angular momentums for the two primitives of 1 and 2
	 * See the document of Gaussian Primitive Product Theorem for more information
	 * \param k   order for angular momentum of x 
	 * \param l1  exponent for PAx
	 * \param l2  exponent for PBx
	 * \param x1  PAx
	 * \param x2  PBx
	 * \return    fk(PAx,PBx,l1,l2)
	 */
	Double Fk(const UInt& k, const UInt& l1, const UInt& l2, 
			const Double& x1, const Double& x2);

	/**
	 * This function is used to calculate the overlap integral 
	 * over two Gaussian primitive functions i and j
	 * ia is the exponential factor for i primitive function
	 * ja is the exponential factor for j primitive function
	 * ic is the coefficient for i primitive function
	 * jc is the coefficient for j primitive function
	 * A  is the center for i primitive function
	 * B  is the center for j primitive function
	 * AB2 is the distance between two primitive centers
	 *
	 * \param d         : exp(-ia*ja/(ia+ja)*AB2)*ic*jc
	 * \param li, mi, ni: angular momentum for i primitive (x^ly^mz^n)
	 * \param lj, mj, nj: angular momentum for j primitive (x^ly^mz^n)
	 * \param alpha     : ia + ja
	 * \param isSameCenter: wether the two primitives share the same center?
	 * \return the overlap integral
	 */
	Double overlapIntegral(const Double& d,   const Double& alpha, 
			const Double& PAx, const Double& PAy, const Double& PAz,  
			const Double& PBx, const Double& PBy, const Double& PBz,  
			const UInt& li, const UInt& mi, const UInt& ni,
			const UInt& lj, const UInt& mj, const UInt& nj, bool isSameCenter);

}

#endif
