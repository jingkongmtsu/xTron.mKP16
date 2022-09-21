/**
 * \file lapack.cpp
 * \author fenglai liu
 *
 * \note
 *
 * right now we have two meigen functions, one is dsyevr(old one), and the new one
 * is dsyevd. It seems that the current dsyevd is not giving correct result for 
 * C atom MO (HF G3LARGE), although the MO gives the same energy. However, it makes
 * the following molecule calculation with atomic guess (CH3CH2O) failed to converge
 * to correct value.
 *
 * We need to see why we have this. Currently I make the dsyevd function call to
 * meigen2 so that to disable it.
 */
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include "excep.h"
#include "blas.h"
#include "blas1.h"
#include "scalablevec.h"
#include "omp_setting.h"
#include "lapack.h"
using namespace std;
using namespace excep;
using namespace omp_setting;
using namespace blas;
using namespace lapack;

void lapack::mllt(Double* M, UInt len, UInt ld, bool withOMP) 
{
#if defined(USE_MKL_LAPACK) || defined(USE_GENERAL_LAPACK)
	Int infor = -1;
	char c = 'L';
	Int len_ = static_cast<Int>(len);
	Int ld_  = static_cast<Int>(ld);

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	spotrf(&c, &len_, M, &ld_, &infor);
#else
	dpotrf(&c, &len_, M, &ld_, &infor);

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

#endif
	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","mllt",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}
#else
	string info = "macro definition is wrong"; 
	Excep excep("lapack","mllt",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,info);
	handleExcep(excep);
#endif
}

void lapack::inverseSPDMatrix(Double* M, UInt len, bool withOMP) 
{
#if defined(USE_MKL_LAPACK) || defined(USE_GENERAL_LAPACK)

	// first step, use the mllt to do factorization
	Int infor = -1;
	char c = 'L';
	Int len_ = static_cast<Int>(len);
	Int ld_  = len_;

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	spotrf(&c, &len_, M, &ld_, &infor);
#else
	dpotrf(&c, &len_, M, &ld_, &infor);
#endif
	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","inverseSPDMatrix::potrf",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

	// secondly, do the inverse
	// reset the infor
	infor = -1;

#ifdef WITH_SINGLE_PRECISION
	spotri(&c, &len_, M, &ld_, &infor);
#else
	dpotri(&c, &len_, M, &ld_, &infor);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","inverseSPDMatrix::potri",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}
#else
	string info = "macro definition is wrong"; 
	Excep excep("lapack","inverseSPDMatrix",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,info);
	handleExcep(excep);
#endif
}

void lapack::inverseGeneralMatrix(Double* M, UInt len, bool withOMP) 
{
#if defined(USE_MKL_LAPACK) || defined(USE_GENERAL_LAPACK)

	// firstly, do the LU decomposition
	Int infor = -1;
	Int len_ = static_cast<Int>(len);
	IntVec ipiv(len_);

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	sgetrf(&len_,&len_,M,&len_,&ipiv.front(),&infor);
#else
	dgetrf(&len_,&len_,M,&len_,&ipiv.front(),&infor);
#endif
	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","inverseMatrix: LU decomposition",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

	// now do the inverse work 
	// check the memory useage first
	Int lwork = -1;
	Double workopt;
#ifdef WITH_SINGLE_PRECISION
	sgetri(&len_,M,&len_,&ipiv.front(),&workopt,&lwork,&infor);
#else
	dgetri(&len_,M,&len_,&ipiv.front(),&workopt,&lwork,&infor);
#endif

	// now do the real work
	lwork = (Int)workopt;
	DoubleVec work(lwork);
	infor = -1;
#ifdef WITH_SINGLE_PRECISION
	sgetri(&len_,M,&len_,&ipiv.front(),&work.front(),&lwork,&infor);
#else
	dgetri(&len_,M,&len_,&ipiv.front(),&work.front(),&lwork,&infor);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","inverseMatrix: getri",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

#else
	string info = "macro definition is wrong"; 
	Excep excep("lapack","inverseMatrix",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,info);
	handleExcep(excep);
#endif
}

void lapack::svd(Double* M, Double* U, Double* V, Double* S, UInt row, 
		UInt col, UInt ld, bool withOMP) 
{
#if defined(USE_MKL_LAPACK) || defined(USE_GENERAL_LAPACK)

	// dimensions
	char job = 'A';
	Int infor= -1;
	Int ldU  = static_cast<Int>(row);
	Int ldV  = static_cast<Int>(col);
	Int row_ = static_cast<Int>(row);
	Int col_ = static_cast<Int>(col);
	Int ld_  = static_cast<Int>(ld);

	// we need to check out how much memory used here
	Int minrowcol = min(row_,col_);
	Double workopt;
	Int lwork = -1;
	IntVec iwork(8*minrowcol);
#ifdef WITH_SINGLE_PRECISION
	sgesdd(&job,&row_,&col_,M,&ld_,S,U,&ldU,V,&ldV,&workopt,
			&lwork,&iwork.front(),&infor);
#else
	dgesdd(&job,&row_,&col_,M,&ld_,S,U,&ldU,V,&ldV,&workopt,
			&lwork,&iwork.front(),&infor);
#endif
	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","svd",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

	// real working code
	lwork = (Int)workopt;
	DoubleVec work(lwork);

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	sgesdd(&job,&row_,&col_,M,&ld_,S,U,&ldU,V,&ldV,&work.front(),
			&lwork,&iwork.front(),&infor);
#else
	dgesdd(&job,&row_,&col_,M,&ld_,S,U,&ldU,V,&ldV,&work.front(),
			&lwork,&iwork.front(),&infor);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","svd",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}
#else
	string info = "macro definition is wrong"; 
	Excep excep("lapack","svd",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,info);
	handleExcep(excep);
#endif
}

void lapack::meigen2(Double* M, Double* V, UInt n, bool withOMP) 
{
#if defined(USE_MKL_LAPACK) || defined(USE_GENERAL_LAPACK)

	// basic parameters
	// since we calculate all of eigenvalues and eigenvectors
	// thus we do not set the range, use null instead
	Int infor = -1;
	char jobz = 'V'; 
	char uplo = 'L';
	Int ld  = static_cast<Int>(n);
	Int n_  = static_cast<Int>(n);

	// working space, prepare for work space query
	Int lwork  = -1;
	Int liwork = -1;
	Double workopt;
	Int iworkopt;

	// memory check up
#ifdef WITH_SINGLE_PRECISION
	ssyevd(&jobz, &uplo, &n_, M, &ld, V, &workopt, &lwork, &iworkopt, &liwork, &infor);
#else
	dsyevd(&jobz, &uplo, &n_, M, &ld, V, &workopt, &lwork, &iworkopt, &liwork, &infor);
#endif
	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","meigen with syevd",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

	// real working code
	lwork = (Int)workopt;
	liwork = (Int)iworkopt;
	DoubleVec work(lwork);
	IntVec iwork(liwork);

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	ssyevd(&jobz, &uplo, &n_, M, &ld, V, &work.front(), &lwork, &iwork.front(), &liwork, &infor);
#else
	dsyevd(&jobz, &uplo, &n_, M, &ld, V, &work.front(), &lwork, &iwork.front(), &liwork, &infor);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","meigen with syevd",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

#else
	string info = "macro definition is wrong"; 
	Excep excep("lapack","meigen with syevd",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,info);
	handleExcep(excep);
#endif
}

void lapack::meigen(Double* M, Double* V, UInt n, bool withOMP) 
{
#if defined(USE_MKL_LAPACK) || defined(USE_GENERAL_LAPACK)

	// basic parameters
	// since we calculate all of eigenvalues and eigenvectors
	// thus we do not set the range, use null instead
	Int infor = -1;
	char jobz = 'V'; 
	char range = 'A';
	char uplo = 'L';
	Int ld = static_cast<Int>(n);
	Double nullD = ZERO;
	Int nullI = 0;
	Double abstol = ZERO;

	// working space, prepare for work space query
	Int lwork = -1;
	Int liwork = -1;
	Double workopt;
	Int iworkopt;

	// output parameters
	Int neigenvalues = 0;
	DoubleVec z(n*n);  // stores the output eigenvectors
	Int ldz = static_cast<Int>(n);
	Int n_  = static_cast<Int>(n);
	IntVec isuppz(2*n); 
#ifdef WITH_SINGLE_PRECISION
	ssyevr(&jobz, &range, &uplo, &n_, M, &ld, &nullD, &nullD, &nullI, &nullI, 
			&abstol, &neigenvalues, V, &z.front(), &ldz, &isuppz.front(), 
			&workopt, &lwork, &iworkopt, &liwork, &infor);
#else
	dsyevr(&jobz, &range, &uplo, &n_, M, &ld, &nullD, &nullD, &nullI, &nullI, 
			&abstol, &neigenvalues, V, &z.front(), &ldz, &isuppz.front(), 
			&workopt, &lwork, &iworkopt, &liwork, &infor);
#endif
	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","meigen with syevr",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

	// real working code
	lwork = (Int)workopt;
	liwork = (Int)iworkopt;
	DoubleVec work(lwork);
	IntVec iwork(liwork);

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	ssyevr(&jobz, &range, &uplo, &n_, M, &ld, &nullD, &nullD, &nullI, &nullI, 
			&abstol, &neigenvalues, V, &z.front(), &ldz, &isuppz.front(), 
			&work.front(), &lwork, &iwork.front(), &liwork, &infor);
#else
	dsyevr(&jobz, &range, &uplo, &n_, M, &ld, &nullD, &nullD, &nullI, &nullI, 
			&abstol, &neigenvalues, V, &z.front(), &ldz, &isuppz.front(), 
			&work.front(), &lwork, &iwork.front(), &liwork, &infor);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","meigen with syevr",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

	// now overwrite the M
	vcopy(&z.front(),M,n*n);

#else
	string info = "macro definition is wrong"; 
	Excep excep("lapack","meigen with syevr",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,info);
	handleExcep(excep);
#endif
}

void lapack::symAXeqB(const char& state, Double* A, UInt n, Double* B, UInt nRHS, bool withOMP)
{
#if defined(USE_MKL_LAPACK) || defined(USE_GENERAL_LAPACK)

#ifdef DEBUG
	if (state != 'U' && state != 'L') {
		string infor = "side is invalid!"; 
		Excep excep("lapack","symAXeqB",EXCEPTION_CHAR_VARIABLE_INVALID,infor);
		handleExcep(excep);
	}
#endif

	// basic parameters
	Int infor = -1;
	char uplo = state;
	Int ldA   = static_cast<Int>(n);
	Int ldB   = static_cast<Int>(n);
	Int n_    = static_cast<Int>(n);
	Int n_RHS = static_cast<Int>(nRHS);
	IntVec ipiv(n,0);

	// working space, prepare for work space query
	Int lwork = -1;
	Double workopt;

	// output parameters
#ifdef WITH_SINGLE_PRECISION
	ssysv(&uplo,&n_,&n_RHS,A,&ldA,&ipiv.front(),B,&ldB,&workopt,&lwork,&infor);
#else
	dsysv(&uplo,&n_,&n_RHS,A,&ldA,&ipiv.front(),B,&ldB,&workopt,&lwork,&infor);
#endif
	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","symAXeqB",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

	// real working code
	lwork = (Int)workopt;
	DoubleVec work(lwork);

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	ssysv(&uplo,&n_,&n_RHS,A,&ldA,&ipiv.front(),B,&ldB,&work.front(),&lwork,&infor);
#else
	dsysv(&uplo,&n_,&n_RHS,A,&ldA,&ipiv.front(),B,&ldB,&work.front(),&lwork,&infor);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","symAXeqB",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

#else
	string info = "macro definition is wrong"; 
	Excep excep("lapack","symAXeqB",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,info);
	handleExcep(excep);
#endif

}

void lapack::AXeqB(Double* A, UInt n, Double* B, UInt nRHS, bool withOMP)
{
#if defined(USE_MKL_LAPACK) || defined(USE_GENERAL_LAPACK)

	// basic parameters
	Int infor = -1;
	Int ldA   = static_cast<Int>(n);
	Int ldB   = static_cast<Int>(n);
	Int n_    = static_cast<Int>(n);
	Int n_RHS = static_cast<Int>(nRHS);
	IntVec ipiv(n,0);

	// set omp to run
	if (withOMP) {
		omp_init();
	}

#ifdef WITH_SINGLE_PRECISION
	sgesv(&n_,&n_RHS,A,&ldA,&ipiv.front(),B,&ldB,&infor);
#else
	dgesv(&n_,&n_RHS,A,&ldA,&ipiv.front(),B,&ldB,&infor);
#endif

	// turn off omp
	if (withOMP) {
		omp_turnoff();
	}

	if (infor != 0) {
		string info = "LAPACK vendor's function returns fail message: " 
			+ boost::lexical_cast<string>(infor); 
		Excep excep("lapack","AXeqB",LAPACK_FUNCTION_FAILED,info);
		handleExcep(excep);
	}

#else
	string info = "macro definition is wrong"; 
	Excep excep("lapack","AXeqB",MACRO_IN_MATH_NOT_PROPERLY_DEFINED,info);
	handleExcep(excep);
#endif

}
