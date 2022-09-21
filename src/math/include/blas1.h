/**
 * \file   blas1.h
 * \author Fenglai Liu
 *
 * the BLAS1 level of functions are separated out from blas.h
 *
 * this is because we want to enforce the possible use of Phi co-processor 
 * on the BLAS1 library
 */
#ifndef BLAS1_H
#define BLAS1_H

// begin the offload attribute region
//#pragma offload_attribute (push, target (mic))

#include "libgen.h"

namespace blas {

	/**
	 * search the maximum value of |x(i)-y(i)|
	 */
	Double maxSearch(const Double* x, const Double* y, UInt len);

	/**
	 * search the maximum value of |x(i)|
	 */
	Double maxSearch(const Double* x, UInt len);

	/**
	 * make the value of x(i) becomes |x(i)|
	 */
	void vabs(Double* x, UInt len);

	/**
	 * compute the root-mean-square deviation for the given error array
	 * we suppose each term already has the x(i)-y(i)
	 * then result is:
	 *
	 * \f E = \sqrt {\frac{\sum_{i}^{n}A_{i}^2}{n}} \f$
	 */
	Double rmsd(const Double* A, UInt len);

	/**
	 * compute the root-mean-square deviation for the given two vectors
	 * then result is:
	 * \f E = \sqrt {\frac{\sum_{i}^{n}(x_{i}-y_{i})^2}{n}} \f$
	 */
	Double rmsd(const Double* x, const Double* y, UInt len);

	/**
	 * sum_{i} x_{i}*y_{i}
	 */
	Double vdot(const Double* x, const Double* y, UInt len); 

	/**
	 * y = a*x + y
	 */
	void vaxpy(const Double* x, Double* y,  Double a, UInt len);

	/**
	 * copy vector x to y
	 */
	void vcopy(const Double* x, Double* y,  UInt len);

	/**
	 * x = a*x
	 */
	void vscal(Double* x,  Double a,  UInt len); 

	/**
	 * set the array to the value of a
	 */
	void vset(Double* x,  Double a, UInt len); 

	/**
	 * z(i) = x(i)*y(i)
	 */
	void vmul(const Double* x,  const Double* y, Double* z,  UInt len); 

	/**
	 * z(i) = x(i)+y(i)
	 */
	void vadd(const Double* x,  const Double* y, Double* z,  UInt len); 

	/**
	 * z(i) = x(i)-y(i)
	 */
	void vsub(const Double* x,  const Double* y, Double* z,  UInt len);

	/**
	 * y(i) = 1/x(i)
	 */
	void vinv(const Double* x, Double* y,  UInt len); 

	/**
	 * z(i) = x(i)/y(i)
	 */
	void vdiv(const Double* x,  const Double* y, Double* z,  UInt len);

	/**
	 * w(i) = x(i) + y(i) + z(i) + w(i)
	 */
	void v3add(const Double* x,  const Double* y,  const Double* z, Double* w,  UInt len); 

	/**
	 * sum_{i} x_{i}*y_{i}*z(i)
	 */
	Double vdot3(const Double* x, const Double* y,  const Double* z,  UInt len);

	/**
	 * sum_{i} w_{i}*x_{i}*y_{i}*z(i)
	 */
	Double vdot4(const Double* w, const Double* x, const Double* y,  const Double* z,  UInt len);

	/**
	 * z(i) = x(i)y(i) + z(i)
	 */
	void vmuladd(const Double* x, const Double* y, Double* z,  UInt len);

	/**
	 * z(i) = a*x(i)y(i) + z(i)
	 */
	void vaxyaddz(const Double* x, const Double* y,  Double a, Double* z,  UInt len);

	/**
	 * z(i) = x(i)y(i) - z(i)
	 */
	void vmulsub(const Double* x,  const Double* y, Double* z,  UInt len);

	/**
	 * w(i) = x(i)*y(i)*z(i) + w(i)
	 */
	void vmul3add(const Double* x,  const Double* y, const Double* z, Double* w, UInt len);

	/**
	 * w(i) = c1*x(i) + c2*y(i) + c3*z(i) + w(i)
	 */
	void v3axpy(const Double& c1, const Double& c2, const Double& c3,
			const Double* x, const Double* y, const Double* z, Double* w, UInt len);

}

// end the offload attribute region
//#pragma offload_attribute (pop)

#endif

