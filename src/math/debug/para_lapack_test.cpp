//
// April 2014
// test the following functions for intel library:
// CD decomposition, SVD, inverse(LU/CD),eigen, linear equation
// only on timing 
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<string> 
#include<vector>
#include "blas.h"
#include "blas1.h"
#include "lapack.h"
#include "tbb/tbb.h"
using namespace blas;
using namespace lapack;
using namespace std;
using namespace tbb;

void para_lapack_test()
{
	// lapcck testing
	cout << "----------------------------------------------" << endl;
	cout << "parallel lapack testing                       " << endl;
	cout << "Cholesky decomposition                        " << endl;
	cout << "SVD decomposition                             " << endl;
	cout << "Symmetrical eigenvectors solver               " << endl;
	cout << "inverse a square matrix                       " << endl;
	cout << "we do not test correctness, only timing will  " << endl;
	cout << "tested here                                   " << endl;
	cout << "----------------------------------------------" << endl;
	bool cdtest = true;
	bool svdtest = true;
	bool eigentest = true;
	bool linearEqSolver = true;
	bool inverse = true;

	// we also test the switch on/off for the multi-threads mode
	for(UInt times=0; times<4; times++) {

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

		// cholesky decomposition
		if (cdtest) {

			cout << endl << endl;
			cout << "----------------------------------------------" << endl;
			cout << "Cholesky decomposition                        " << endl;
			cout << "----------------------------------------------" << endl;
			// basic setting for testing matrix
			UInt rowa = 10000;
			UInt cola = 10000;

			// set up matrix data
			// we need to set up a full rank matrix
			vector<Double> a(rowa*cola,ONE);
			for(UInt i=0; i<rowa; i++) {
				a[i+rowa*i] = rowa;
			}

			// now in comparison
			tick_count t0 = tick_count::now();
			mllt(&a.front(),rowa,rowa,doPara);
			tick_count t1 = tick_count::now();
			Double t = (t1-t0).seconds();
			printf("%s  %-12.6f\n", "mllt(10000*10000) time with seconds ", t);
		}

		//testing SVD
		if (svdtest) {

			cout << endl << endl;
			cout << "----------------------------------------------" << endl;
			cout << "SVD                                           " << endl;
			cout << "----------------------------------------------" << endl;

			// define the data testing matrix
			UInt row = 5000;
			UInt col = 5000;
			vector<Double> x(row*col,ONE);
			for(UInt j=0; j<col; j++) {
				x[j+j*row] = row;
			}

			vector<Double> S(row*row);
			vector<Double> VT(col*col);
			vector<Double> eigen(row);
			tick_count t0 = tick_count::now();
			svd(&x.front(),&S.front(),&VT.front(),&eigen.front(),row,col,row,doPara);
			tick_count t1 = tick_count::now();
			Double t = (t1-t0).seconds();
			printf("%s  %-12.6f\n", "svd(5000*5000) time with seconds ", t);
		}

		// eigen test
		if (eigentest) {

			cout << endl << endl;
			cout << "----------------------------------------------" << endl;
			cout << "EIGEN decomposition                           " << endl;
			cout << "----------------------------------------------" << endl;
			// set up the original matrix
			UInt n = 10000;
			vector<Double> b(n*n,ONE);
			for(UInt i=0; i<n; i++) {
				b[i+i*n] = n; 
			}

			// real work
			vector<Double> e(n);
			tick_count t0 = tick_count::now();
			meigen(&b.front(),&e.front(),n,doPara);
			tick_count t1 = tick_count::now();
			Double t = (t1-t0).seconds();
			printf("%s  %-12.6f\n", "meigen(10000*10000) time with seconds ", t);
		}

		// test linear equation solver
		if(linearEqSolver) {

			cout << endl << endl;
			cout << "----------------------------------------------" << endl;
			cout << "AX = B solver                                 " << endl;
			cout << "----------------------------------------------" << endl;
			UInt n = 10000;
			vector<Double> A(n*n,ONE);
			for(UInt i=0; i<n; i++) {
				A[i+n*i] = n;
			}
			UInt rowB = n;
			UInt colB = 10;
			vector<Double> B(rowB*colB,0.5E0);

			tick_count t0 = tick_count::now();
			symAXeqB('L', &A.front(), rowB, &B.front(), colB,doPara);
			tick_count t1 = tick_count::now();
			Double t = (t1-t0).seconds();
			printf("%s  %-12.6f\n", "symAXeqB time with seconds ", t);
			printf("AX=B, A is in dimension of (10000*10000), B is (10000*10)\n");

			// now for AX=B assymetrical case
			for(UInt i=0; i<n; i++) {
				A[i+n*i] = n;
			}
			B.assign(rowB*colB,0.5E0);

			tick_count t2 = tick_count::now();
			AXeqB(&A.front(),rowB,&B.front(),colB,doPara);
			tick_count t3 = tick_count::now();
			t = (t3-t2).seconds();
			printf("%s  %-12.6f\n", "AXeqB time with seconds ", t);
			printf("AX=B, A is in dimension of (10000*10000), B is (10000*10)\n");
		}

		if (inverse) {

			cout << endl << endl;
			cout << "----------------------------------------------" << endl;
			cout << "INVERSE MATRIX with LU decomposition          " << endl;
			cout << "----------------------------------------------" << endl;
			// input matrix
			UInt n = 10000;
			vector<Double> A(n*n,ONE);
			for(UInt i=0; i<n; i++) {
				A[i+n*i] = n;
			}

			// now do inverse
			tick_count t0 = tick_count::now();
			inverseGeneralMatrix(&A.front(),n,doPara);
			tick_count t1 = tick_count::now();
			Double t = (t1-t0).seconds();
			printf("%s  %-12.6f\n", "inverse matrix(LU) time with seconds(10000*10000) ", t);

		}

		if (inverse) {

			cout << endl << endl;
			cout << "----------------------------------------------" << endl;
			cout << "INVERSE MATRIX with CD decomposition          " << endl;
			cout << "----------------------------------------------" << endl;
			// input matrix
			UInt n = 10000;
			vector<Double> A(n*n,ONE);
			for(UInt i=0; i<n; i++) {
				A[i+n*i] = n;
			}

			// now do inverse
			tick_count t0 = tick_count::now();
			inverseSPDMatrix(&A.front(),n,doPara);
			tick_count t1 = tick_count::now();
			Double t = (t1-t0).seconds();
			printf("%s  %-12.6f\n", "inverse matrix(CD) time with seconds(10000*10000) ", t);
		}
	}
}
