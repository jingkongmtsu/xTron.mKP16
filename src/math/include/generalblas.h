/**
 * head files when we use general blas from netlib etc.
 * \author Fenglai Liu
 *
 */
#ifndef GENERALBLAS_H
#define GENERALBLAS_H
#include "libgen.h"

// names conversion for double type
#define  dgemv    dgemv_
#define  dsymv    dsymv_
#define  dgemm    dgemm_
#define  dsymm    dsymm_

// names conversion for single floating type
#define  sgemv    sgemv_
#define  ssymv    ssymv_
#define  sgemm    sgemm_
#define  ssymm    ssymm_

extern "C" {

#ifdef WITH_SINGLE_PRECISION

	extern void sgemv(char* symA, Int* rowA, Int* colA, 
			Double* alpha, const Double* A, Int* ldA, const Double* x, 
			Int* incx, Double* beta, Double* y, Int* incy);
	extern void ssymv(char* symA, Int* n, Double* alpha, const Double* A, 
			Int* ldA, const Double* x, Int* incx, Double* beta, Double* y, Int* incy);
	extern void sgemm(char* symA, char* symB, Int* rowa, Int* colb, 
			Int* cola, Double* alpha, const Double* A, Int* ldA, 
			const Double* B, Int* ldB, Double* beta, Double* C, Int* ldC);
	extern void ssymm(char* side, char* uplo, Int* m, Int* n, Double* alpha, 
			const Double* A, Int* ldA, const Double* B, Int* ldB, Double* beta, 
			Double* C, Int* ldC);

#else

	extern void dgemv(char* symA, Int* rowA, Int* colA, 
			Double* alpha, const Double* A, Int* ldA, const Double* x, 
			Int* incx, Double* beta, Double* y, Int* incy);
	extern void dsymv(char* symA, Int* n, Double* alpha, const Double* A, 
			Int* ldA, const Double* x, Int* incx, Double* beta, Double* y, Int* incy);
	extern void dgemm(char* symA, char* symB, Int* rowa, Int* colb, 
			Int* cola, Double* alpha, const Double* A, Int* ldA, 
			const Double* B, Int* ldB, Double* beta, Double* C, Int* ldC);
	extern void dsymm(char* side, char* uplo, Int* m, Int* n, Double* alpha, 
			const Double* A, Int* ldA, const Double* B, Int* ldB, Double* beta, 
			Double* C, Int* ldC);

#endif
}

#endif
