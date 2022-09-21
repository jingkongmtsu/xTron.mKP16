/**
 * \file   lapack.h
 * \author Fenglai Liu 
 *
 */
#ifndef GENERALLAPACK_H
#define GENERALLAPACK_H
#include "libgen.h"

// names conversion
#define  dpotrf  dpotrf_
#define  dpotri  dpotri_
#define  dgesdd  dgesdd_
#define  dsyevr  dsyevr_
#define  dsyevd  dsyevd_
#define  dsysv   dsysv_
#define  dgetrf  dgetrf_
#define  dgetri  dgetri_
#define  dgesv   dgesv_

#define  spotrf  spotrf_
#define  spotri  spotri_
#define  sgesdd  sgesdd_
#define  ssyevr  ssyevr_
#define  ssyevd  ssyevd_
#define  ssysv   ssysv_
#define  sgetrf  sgetrf_
#define  sgetri  sgetri_
#define  sgesv   sgesv_

extern "C" {

#ifdef WITH_SINGLE_PRECISION

	extern void sgetrf(Int* row, Int* col, Double* M, Int* lda, Int* ipiv, Int* infor);
	extern void sgetri(Int* N, Double* A, Int* lda, Int* IPIV, Double* WORK, 
			Int* lwork, Int* INFO);
	extern void spotrf(char* c, Int* len, Double* M, Int* ld, Int* infor);
	extern void spotri(char* c, Int* len, Double* M, Int* ld, Int* infor);
	extern void sgesdd(char* job, Int* row, Int* col,Double* M,
			Int* ld, Double* S, Double* U, Int* ldU, Double* V, Int* ldV,
			Double* work, Int* lwork, Int* iwork, Int* infor);
	extern void ssyevr(char* jobz, char* range, char* uplo, 
			Int* n, Double* M, Int* ld, Double* nullD1, Double* nullD2, 
			Int* nullI1, Int* nullI2, Double* abstol, Int* neigenvalues, 
			Double* V, Double* z, Int* ldz, Int* isuppz, 
			Double* work, Int* lwork, Int* iwork, Int* liwork, Int* infor);
	extern void ssyevd(char* jobz, char* uplo, Int* n, Double* a, Int* lda, Double* w, Double* work, 
			Int* lwork, Int* iwork, Int* liwork, Int* info);
	extern void ssysv(char* uplo, Int* n, Int* nrhs, Double* a, Int* lda, 
			Int* ipiv, Double* b, Int* ldb, Double* work, Int* lwork, Int* info);
	extern void sgesv(Int* n, Int* nrhs, Double* a, Int* lda, 
			Int* ipiv, Double* b, Int* ldb, Int* info);
#else

	extern void dgetrf(Int* row, Int* col, Double* M, Int* lda, Int* ipiv, Int* infor);
	extern void dgetri(Int* N, Double* A, Int* lda, Int* IPIV, Double* WORK, 
			Int* lwork, Int* INFO);
	extern void dpotrf(char* c, Int* len, Double* M, Int* ld, Int* infor);
	extern void dpotri(char* c, Int* len, Double* M, Int* ld, Int* infor);
	extern void dgesdd(char* job, Int* row, Int* col,Double* M,
			Int* ld, Double* S, Double* U, Int* ldU, Double* V, Int* ldV,
			Double* work, Int* lwork, Int* iwork, Int* infor);
	extern void dsyevr(char* jobz, char* range, char* uplo, 
			Int* n, Double* M, Int* ld, Double* nullD1, Double* nullD2, 
			Int* nullI1, Int* nullI2, Double* abstol, Int* neigenvalues, 
			Double* V, Double* z, Int* ldz, Int* isuppz, 
			Double* work, Int* lwork, Int* iwork, Int* liwork, Int* infor);
	extern void dsyevd(char* jobz, char* uplo, Int* n, Double* a, Int* lda, Double* w, Double* work, 
			Int* lwork, Int* iwork, Int* liwork, Int* info);
	extern void dsysv(char* uplo, Int* n, Int* nrhs, Double* a, Int* lda, 
			Int* ipiv, Double* b, Int* ldb, Double* work, Int* lwork, Int* info);
	extern void dgesv(Int* n, Int* nrhs, Double* a, Int* lda, 
			Int* ipiv, Double* b, Int* ldb, Int* info);
#endif
}

#endif
