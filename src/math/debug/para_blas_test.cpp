//
// note:
// April 2014:
// test vdot and vaxpy
// April 2014
// test BLAS2 and BLAS3 for intel library
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<string> 
#include<vector>
#include<cmath>
#include "libgen.h"
#include "tbb/tbb.h"
#include "blas.h"
using namespace blas;
using namespace tbb;
using namespace std;

void para_blas_test(bool testingBLAS1, bool testingBLAS2, bool testingBLAS3)
{
	// blas2 testing
	if (testingBLAS2) {

		cout << "----------------------------------------------" << endl;
		cout << "parallel BLAS2 testing with openmp            " << endl;
		cout << "all matrix are in column major                " << endl;
		cout << "we do not test the correctness of the result  " << endl;
		cout << "since the parallel version is totally provided" << endl;
		cout << "by the BLAS vendor, therefore the correctness " << endl;
		cout << "of parallel code will be guaranteed by them   " << endl;
		cout << "----------------------------------------------" << endl;

		// basic setting for testing matrix
		UInt rowa = 50000;
		UInt cola = 50000;

		// set up matrix and vector
		// we need to set up a full rank matrix
		vector<Double> a(rowa*cola);

		// set up two vector
		vector<Double> x(rowa);
		vector<Double> y(rowa);

		for(UInt times=0; times<6; times++) {

			// set the parallel choice
			bool doPara = false;
			if (times%2 == 1) doPara = true;
			cout << "********************************" << endl;
			if (doPara) {
				cout << "do para   section               " << endl;               
			}else{
				cout << "do serial section               " << endl;               
			}
			cout << "********************************" << endl;

			// symmetrical testing
			a.assign(a.size(), ONE);
			x.assign(rowa, ONE);
			y.assign(rowa, ZERO);
			tick_count t0 = tick_count::now();
			vmsymv(&a.front(), rowa, rowa, &x.front(), &y.front(), ONE, ONE, doPara);
			tick_count t1 = tick_count::now();
			Double t = (t1-t0).seconds();
			printf("%s  %-12.6f\n", "vmsymv time with seconds ", t);

			// general testing
			a.assign(a.size(), ONE);
			x.assign(rowa, ONE);
			y.assign(rowa, ZERO);
			tick_count t2 = tick_count::now();
			vmgemv(&a.front(), 'N', rowa, rowa, rowa, &x.front(), &y.front(), ONE, ONE, doPara);
			tick_count t3 = tick_count::now();
			t = (t3-t2).seconds();
			printf("%s  %-12.6f\n", "vmgemv time with seconds ", t);
		}
	}

	// blas3 testing
	if (testingBLAS3) {

		cout << "----------------------------------------------" << endl;
		cout << "parallel BLAS3 testing with openmp            " << endl;
		cout << "all matrix are in column major                " << endl;
		cout << "we do not test the correctness of the result  " << endl;
		cout << "since the parallel version is totally provided" << endl;
		cout << "by the BLAS vendor, therefore the correctness " << endl;
		cout << "of parallel code will be guaranteed by them   " << endl;
		cout << "----------------------------------------------" << endl;

		// basic setting for testing matrix
		UInt rowa = 5000;
		UInt cola = 5500;

		// set up matrix data
		vector<Double> A(rowa*cola,ONE);

		// first, B is symmetrical and multiply on the left side
		UInt rowb = rowa;
		UInt colb = rowa;
		vector<Double> B(rowb*colb,ONE);
		UInt rowc = rowb;
		UInt colc = cola;
		vector<Double> C(rowc*colc);

		for(UInt times=0; times<4; times++) {

			// set the parallel choice
			bool doPara = false;
			if (times%2 == 0) doPara = true;
			cout << "********************************" << endl;
			if (doPara) {
				cout << "do para   section               " << endl;               
			}else{
				cout << "do serial section               " << endl;               
			}
			cout << "********************************" << endl;

			tick_count t0 = tick_count::now();
			msymmul(&B.front(),&A.front(),&C.front(),rowb,cola,'L', ONE, ZERO, doPara);
			tick_count t1 = tick_count::now();
			Double t = (t1-t0).seconds();
			cout << "testing L side symmetrical matrix multiplication: C = B*A " << endl;
			printf("%s  %-12.6f\n", "msymmul time with seconds ", t);

			//
			// now test the mmul, without LD version
			// however, we will try to do it with transpose
			//
			C.assign(C.size(),ONE);
			tick_count t2 = tick_count::now();
			mmul(&A.front(),&B.front(),&C.front(),rowa,cola,rowb,colb,'T','N',ONE,ONE,doPara);
			tick_count t3 = tick_count::now();
			t = (t3-t2).seconds();
			cout << "testing matrix multiplication: C = A*B^{T}+C " << endl;
			printf("%s  %-12.6f\n", "mmul time with seconds ", t);
		}
	}
}
