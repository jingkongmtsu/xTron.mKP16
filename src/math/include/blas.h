/**
 * \file   blas.h
 * \author Fenglai Liu
 *
 * We provide BLAS 2/3 level of operations
 *
 * the blas1 level of code is contained in blas1.h
 */
#ifndef BLAS_H
#define BLAS_H
#include "libgen.h"

#if defined(USE_MKL)
#define USE_MKL_BLAS
#include "mkl.h"
#elif defined(USE_GENERAL_BLAS)
#include "generalblas.h"
#endif

namespace blas {

	/*********************************************************************************
	 *                               BLAS2 LIBRARY                                   *
	 *********************************************************************************/
	/**
	 * y = \alpha op(A) x + \beta y
	 * A is the general matrix, x and y are vectors
	 * this is the wrapper function for dgemv/sgemv
	 */
	void vmgemv(const Double* A, char symA,  UInt rowA,  UInt colA,
			UInt ldA, const Double* x, Double* y,  Double alpha,  Double beta,  
			bool withOMP = false);

	/**
	 * y = \alpha Ax + \beta y
	 * A is the symmetrical matrix in full storage, x and y are vectors
	 * this is the wrapper function for dsymv/ssymv
	 */
	void vmsymv(const Double* A,  UInt n, UInt ldA, const Double* x, 
			Double* y,  Double alpha,  Double beta, bool withOMP = false);

	/*********************************************************************************
	 *                               BLAS3 LIBRARY                                   *
	 *             we note that here it does not support pack form matrix            *
	 *********************************************************************************/

	/**
	 * matrix multiplication
	 * C = alpha*op(A)*op(B) + beta*C
	 * op(A) or op(B) mean A could be normal or transpose matrix
	 * alpha and beta are two numbers
	 * this is the complete version of matrix multiplication
	 * The transpose wll be taken care in the code, in the interface we just pass
	 * the normal dimension without taking care of the transpose
	 *
	 * \param   ldA: leading dimension of A
	 * \param   ldB: leading dimension of B
	 * \param   ldC: leading dimension of C
	 * \param  symA: symbol to show the state of A, could be only "N" or "T"(transpose)
	 * \param  symB: symbol to show the state of B, could be only "N" or "T"(transpose)
	 *
	 */
	void mmul(const Double* A, const Double* B, Double* C,  UInt rowA,
			UInt colA,  UInt rowB,  UInt colB, UInt ldA,  UInt ldB,  UInt ldC,
			char symA, char symB,  Double alpha,  Double beta, bool withOMP = false); 

	/**
	 * if the matrix is not the submatrix, then the leading dimension is just the row
	 * number. So this short version is taking advantage of this, we do not
	 * have the leading dimension anymore. However, we should take care fo ldC. 
	 * For LdC, it depends on A is transpose or not
	 */
	void mmul(const Double* A, const Double* B, Double* C,  UInt rowA,
			UInt colA,  UInt rowB,  UInt colB, char symA, char symB, 
			Double alpha, Double beta, bool withOMP = false); 

	/**
	 * furthermore, if both of A and B are all in 'N' state, then we can even drop 
	 * the symA and symB
	 */
	void mmul(const Double* A, const Double* B, Double* C,  UInt rowA,
			UInt colB,  UInt colA, Double alpha, Double beta, bool withOMP = false);

	/**
	 * matrix multiplication
	 * C = alpha*A*B + beta*C or C = alpha*B*A + beta*C
	 * A is symmetrical matrix in dimension N*N, B is M*N or N*M matrix
	 * alpha and beta are two numbers
	 *
	 * \param   ldA: leading dimension of A
	 * \param   ldB: leading dimension of B
	 * \param   ldC: leading dimension of C
	 * \param  side: A is on the left "L", or on the right "R"
	 *
	 */
	void msymmul(const Double* A, const Double* B, Double* C,  UInt N,
			UInt M,  UInt ldA,  UInt ldB,  UInt ldC, char side,  Double alpha,  
			Double beta, bool withOMP = false);

	/**
	 * This is the msymmul without leading dimension
	 *
	 */
	void msymmul(const Double* A, const Double* B, Double* C,  UInt N,
			UInt M, char side,  Double alpha,  Double beta, bool withOMP = false);

}

#endif

