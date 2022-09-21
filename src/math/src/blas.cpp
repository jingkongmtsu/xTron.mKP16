/**
 * \file blas.cpp
 * \author fenglai liu 
 *
 * functions for BLAS 2 and 3
 * For BLAS1, the code is in blas1.cpp
 *
 * for BLAS2 and BLAS3, we need to call the vendor's function to do the job
 * the parallization is responsible by blas vendor
 *
 * ALL of functions in BLAS1 and BLAS2 has increment of 1 (see blas2, incx and incy)
 */
#include <cstdio>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include "excep.h"
#include "blas.h"
#include "omp_setting.h"
using namespace excep;
using namespace omp_setting;
using namespace blas;

/*********************************************************************************
 *                     ####      BLAS2 LIBRARY                                   *
 *********************************************************************************/
void blas::vmgemv(const Double* A, char symA,  UInt rowA,  UInt colA,
			 UInt ldA,  const Double* x, Double* y,  Double alpha,  Double beta, bool withOMP) 
{
	Int row_A  = static_cast<Int>(rowA);
	Int col_A  = static_cast<Int>(colA);
	Int ld_A   = static_cast<Int>(ldA);
	Int inc_x  = 1;
	Int inc_y  = 1;

#if defined(USE_MKL_BLAS) || defined(USE_GENERAL_BLAS)

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	sgemv(&symA, &row_A, &col_A, &alpha, A, &ld_A, x, &inc_x, &beta, y, &inc_y);
#else
	dgemv(&symA, &row_A, &col_A, &alpha, A, &ld_A, x, &inc_x, &beta, y, &inc_y);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

#else

	string infor = "macro definition is wrong"; 
	Excep excep("blas","vmgemv",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,infor);
	handleExcep(excep);

#endif
}

void blas::vmsymv(const Double* A,  UInt n, UInt ldA,  const Double* x, 
		Double* y,  Double alpha,  Double beta, bool withOMP) 
{
	char uplo  = 'L';
	Int n_f    = static_cast<Int>(n);
	Int ld_A   = static_cast<Int>(ldA);
	Int inc_x  = 1;
	Int inc_y  = 1;
#if defined(USE_MKL_BLAS) || defined(USE_GENERAL_BLAS)

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	ssymv(&uplo, &n_f, &alpha, A, &ld_A, x, &inc_x, &beta, y, &inc_y);
#else
	dsymv(&uplo, &n_f, &alpha, A, &ld_A, x, &inc_x, &beta, y, &inc_y);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

#else
	string infor = "macro definition is wrong"; 
	Excep excep("blas","vmsymv",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,infor);
	handleExcep(excep);
#endif
}

/*********************************************************************************
 *                     ####      BLAS3 LIBRARY                                   *
 *********************************************************************************/
void blas::mmul(const Double* A, const Double* B, Double* C,  UInt rowA,
			 UInt colA,  UInt rowB,  UInt colB, UInt ldA,  UInt ldB,  UInt ldC,
			char symA, char symB,  Double alpha,  Double beta, bool withOMP) 
{

	// check the input char
#ifdef DEBUG
	if (symA != 'N' && symA != 'T') {
		string infor = "symA is invalid!"; 
		Excep excep("blas","mmul",EXCEPTION_CHAR_VARIABLE_INVALID,infor);
		handleExcep(excep);
	}
	if (symB != 'N' && symB != 'T') {
		string infor = "symB is invalid!"; 
		Excep excep("blas","mmul",EXCEPTION_CHAR_VARIABLE_INVALID,infor);
		handleExcep(excep);
	}
#endif

	// dimension change and check
	UInt rowa = rowA;
	UInt cola = colA;
	UInt colb = colB;
	if (symA == 'T') {
		rowa = colA;
		cola = rowA;
	}
	if (symB == 'T'){
		colb = rowB;
	}
#ifdef DEBUG
	UInt rowb = rowB;
	if (symB == 'T'){
		rowb = colB;
	}
	if(cola != rowb){
		string infor = "row A and col B do not match";
		Excep excep("blas","mmul",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif

	// now call vendor's function
	Int row_a  = static_cast<Int>(rowa);
	Int col_a  = static_cast<Int>(cola);
	Int col_b  = static_cast<Int>(colb);
	Int ld_A   = static_cast<Int>(ldA);
	Int ld_B   = static_cast<Int>(ldB);
	Int ld_C   = static_cast<Int>(ldC);
#if defined(USE_MKL_BLAS) || defined(USE_GENERAL_BLAS)

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	sgemm(&symA, &symB, &row_a, &col_b, &col_a, &alpha, A, &ld_A, B, &ld_B, &beta, C, &ld_C);
#else
	dgemm(&symA, &symB, &row_a, &col_b, &col_a, &alpha, A, &ld_A, B, &ld_B, &beta, C, &ld_C);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

#else
	string infor = "macro definition is wrong"; 
	Excep excep("blas","mmul",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,infor);
	handleExcep(excep);
#endif
}	

void blas::mmul(const Double* A, const Double* B, Double* C,  UInt rowA,
		UInt colA,  UInt rowB,  UInt colB, char symA, char symB, 
		Double alpha, Double beta, bool withOMP) 
{

	// check input char
#ifdef DEBUG
	if (symA != 'N' && symA != 'T') {
		string infor = "symA is invalid!"; 
		Excep excep("blas","mmul",EXCEPTION_CHAR_VARIABLE_INVALID,infor);
		handleExcep(excep);
	}
	if (symB != 'N' && symB != 'T') {
		string infor = "symB is invalid!"; 
		Excep excep("blas","mmul",EXCEPTION_CHAR_VARIABLE_INVALID,infor);
		handleExcep(excep);
	}
#endif

	// dimension handling
	UInt ldA  = rowA;
	UInt ldB  = rowB;
	UInt rowa = rowA;
	UInt cola = colA;
	UInt colb = colB;
	if (symA == 'T') {
		rowa = colA;
		cola = rowA;
	}
	if (symB == 'T'){
		colb = rowB;
	}
	UInt ldC  = rowa;
#ifdef DEBUG
	UInt rowb = rowB;
	if (symB == 'T'){
		rowb = colB;
	}
	if(cola != rowb){
		string infor = "row A and col B do not match";
		Excep excep("blas","mmul",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif

	// call vendor's function
	Int row_a  = static_cast<Int>(rowa);
	Int col_a  = static_cast<Int>(cola);
	Int col_b  = static_cast<Int>(colb);
	Int ld_A   = static_cast<Int>(ldA);
	Int ld_B   = static_cast<Int>(ldB);
	Int ld_C   = static_cast<Int>(ldC);
#if defined(USE_MKL_BLAS) || defined(USE_GENERAL_BLAS)

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	sgemm(&symA, &symB, &row_a, &col_b, &col_a, &alpha, A, &ld_A, B, &ld_B, &beta, C, &ld_C);
#else
	dgemm(&symA, &symB, &row_a, &col_b, &col_a, &alpha, A, &ld_A, B, &ld_B, &beta, C, &ld_C);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

#else
	string infor = "macro definition is wrong"; 
	Excep excep("blas","mmul",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,infor);
	handleExcep(excep);
#endif
}	

void blas::mmul(const Double* A, const Double* B, Double* C,  UInt rowA,
			 UInt colB,  UInt colA, Double alpha, Double beta, bool withOMP) 
{
	char symA  = 'N';
	char symB  = 'N';	   	
	Int row_a  = static_cast<Int>(rowA);
	Int col_a  = static_cast<Int>(colA);
	Int col_b  = static_cast<Int>(colB);
	Int ld_A   = row_a;
	Int ld_B   = col_a;
	Int ld_C   = row_a;
#if defined(USE_MKL_BLAS) || defined(USE_GENERAL_BLAS)

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	sgemm(&symA, &symB, &row_a, &col_b, &col_a, &alpha, A, &ld_A, B, &ld_B, &beta, C, &ld_C);
#else
	dgemm(&symA, &symB, &row_a, &col_b, &col_a, &alpha, A, &ld_A, B, &ld_B, &beta, C, &ld_C);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

#else
	string infor = "macro definition is wrong"; 
	Excep excep("blas","mmul",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,infor);
	handleExcep(excep);
#endif
}	

void blas::msymmul(const Double* A, const Double* B, Double* C,  UInt N,
			 UInt M,  UInt ldA,  UInt ldB,  UInt ldC,
			char side,  Double alpha,  Double beta, bool withOMP) 
{
	// check the input char
#ifdef DEBUG
	if (side != 'R' && side != 'L') {
		string infor = "side is invalid!"; 
		Excep excep("blas","msymmul",EXCEPTION_CHAR_VARIABLE_INVALID,infor);
		handleExcep(excep);
	}
#endif

	// check the dimension
	UInt m = N;
	UInt n = M;
	if (side == 'R'){
		m = M;
		n = N;
	}
	char uplo = 'L';

	// call vendor's function
	Int m_f    = static_cast<Int>(m);
	Int n_f    = static_cast<Int>(n);
	Int ld_A   = static_cast<Int>(ldA);
	Int ld_B   = static_cast<Int>(ldB);
	Int ld_C   = static_cast<Int>(ldC);
#if defined(USE_MKL_BLAS) || defined(USE_GENERAL_BLAS)

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	ssymm(&side, &uplo, &m_f, &n_f, &alpha, A, &ld_A, B, &ld_B, &beta, C, &ld_C);
#else
	dsymm(&side, &uplo, &m_f, &n_f, &alpha, A, &ld_A, B, &ld_B, &beta, C, &ld_C);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

#else
	string infor = "macro definition is wrong"; 
	Excep excep("blas","msymmul",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,infor);
	handleExcep(excep);
#endif
}	

void blas::msymmul(const Double* A, const Double* B, Double* C,  UInt N,
		UInt M, char side,  Double alpha,  Double beta, bool withOMP) 
{
	// check the input char
#ifdef DEBUG
	if (side != 'R' && side != 'L') {
		string infor = "side is invalid!"; 
		Excep excep("blas","msymmul",EXCEPTION_CHAR_VARIABLE_INVALID,infor);
		handleExcep(excep);
	}
#endif

	// check dimension
	UInt m = N;
	UInt n = M;
	if (side == 'R'){
		m = M;
		n = N;
	}
	UInt ldA = N;
	UInt ldB = m;
	UInt ldC = m;
	char uplo = 'L';

	// call vendor's function
	Int m_f    = static_cast<Int>(m);
	Int n_f    = static_cast<Int>(n);
	Int ld_A   = static_cast<Int>(ldA);
	Int ld_B   = static_cast<Int>(ldB);
	Int ld_C   = static_cast<Int>(ldC);
#if defined(USE_MKL_BLAS) || defined(USE_GENERAL_BLAS)

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	ssymm(&side, &uplo, &m_f, &n_f, &alpha, A, &ld_A, B, &ld_B, &beta, C, &ld_C);
#else
	dsymm(&side, &uplo, &m_f, &n_f, &alpha, A, &ld_A, B, &ld_B, &beta, C, &ld_C);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

#else
	string infor = "macro definition is wrong"; 
	Excep excep("blas","msymmul",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,infor);
	handleExcep(excep);
#endif
}	

