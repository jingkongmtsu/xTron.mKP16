/**
 * \file   lapack.h
 * \author Fenglai Liu
 *
 * build the wrapper function to call vendor's lapack functions
 */
#ifndef LAPACK_H
#define LAPACK_H
#include "libgen.h"
#ifdef USE_MKL
#define USE_MKL_LAPACK
#include "mkl.h"
#elif defined(USE_GENERAL_LAPACK)
#include "generallapack.h"
#endif

namespace lapack {

	/**
	 * Computes the Cholesky decomposition of a symmetric positive-definite matrix.
	 * only the lower part of data is used
	 * M = L*L^{T},  L is the lower triangular matrix
	 * we note that L is overwritten into the input M
	 */
	void mllt(Double* M, UInt len, UInt ld, bool withOMP = false); 

	/**
	 * generally inverse a matrix through LU decompostion
	 * this is true for any kind of square matrix
	 * the only thing is that the input matrix should not be 
	 * sigular, else the result matrix is not accurate
	 *
	 * if you use this function to compute inverse, you should
	 * test the result matrix to see whether it's accurate
	 */
	void inverseGeneralMatrix(Double* M, UInt len, bool withOMP = false);

	/**
	 * for a SPD matrix (Symmetric Positive Definite matrix),
	 * it's inverse is much more easier. We can use the 
	 * Cholesky decomposition to do it
	 */
	void inverseSPDMatrix(Double* M, UInt len, bool withOMP = false);

	/**
	 * perform Singular Value Decomposition for a given matrix in m*n
	 * by using the divide and conquor method. 
	 * M = U*S*V^{T}
	 * S is one dimensional array, with length of (min(row,col))
	 * here we note that u should be m*m and it's ld = m
	 * here we note that v should be n*n and it's ld = n
	 * the returned V is not the real V, but it's transpose form V^{T}
	 */
	void svd(Double* M, Double* U, Double* V, Double* S, 
			UInt row, UInt col, UInt ld, bool withOMP = false); 

	/**
	 * Eigen values and eigen vectors under calling syevd function
	 */
	void meigen2(Double* M, Double* V, UInt n, bool withOMP = false);

	/**
	 * Eigen values and eigen vectors under RRR 
	 * algorithm(Relatively Robust Representations), this is recommended by MKL
	 * M is the given symmetric matrix in full storage, and it's overwritten 
	 * by the eigenvectors; the original content only lower part needed
	 * We note that for matrix M, it must be full storage since finally the 
	 * orthogonal eigenvectors would be stored inside
	 * V is the corresponding eigen values.
	 * all of the eigen values should be arranged in ascending way
	 */
	void meigen(Double* M, Double* V, UInt n, bool withOMP = false);

	/**
	 * Computes the solution to the system of linear equations with a real 
	 * symmetric matrix A and multiple right-hand sides B. A is nA*nA
	 * with only upper or lower triangular data(characterized by state), 
	 * B is nA*nRHS. The whole process is to solve:
	 * AX = B 
	 * X is what we are solving. This function actually calls dsysv/ssysv for a
	 * real job.
	 *
	 * Both A and B would be overwritting. B is overwritting with X, 
	 * A is overwritting with is overwritten by the block-diagonal matrix 
	 * D and the multipliers used to obtain the factor U (or L) from 
	 * the factorization of A.
	 *
	 * We note that A should not be singular. Else we will raise an error
	 */
	void symAXeqB(const char& state, Double* A, UInt nA, Double* B, 
			UInt nRHS, bool withOMP = false);

	/**
	 * Computes the solution to the system of linear equations with a real 
	 * matrix A and multiple right-hand sides B. A is nA*nA
	 * B is nA*nRHS. The whole process is to solve:
	 * AX = B 
	 * X is what we are solving. This function actually calls dgesv/sgesv for a
	 * real job.
	 *
	 * Both A and B would be overwritting. B is overwritting with X, 
	 * A is overwritten by the factors L and U from the factorization of 
	 * A = P*L*U; the unit diagonal elements of L are not stored.
	 *
	 */
	void AXeqB(Double* A, UInt nA, Double* B,  UInt nRHS, bool withOMP = false);
}

#endif

