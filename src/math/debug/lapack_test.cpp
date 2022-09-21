//
// for single/double precision (March, 2014, netlib ):
// ALL test are done: CD decomposition, SVD, inverse(LU/CD),eigen, linear equation
//
// April 2014(mkl 11.1 in multi-threading version):
// for LAPACK from intel all functions are tested, single and double
//
// May 2014:
// testing the AX=B for assymetrical case
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<string> 
#include<vector>
#include "blas.h"
#include "lapack.h"
#include "omp_setting.h"
using namespace blas;
using namespace lapack;
using namespace std;
using namespace omp_setting;

extern void lowToUpper(UInt n, vector<Double>& A);

void lapack_test()
{
	// firstly, let's concentrate on serial mode
	// turn openmp off if the BLAS lib use openmp
	omp_turnoff();

	// lapcck testing
	cout << "----------------------------------------------" << endl;
	cout << "lapack testing                                " << endl;
	cout << "Cholesky decomposition                        " << endl;
	cout << "SVD decomposition                             " << endl;
	cout << "Symmetrical eigenvectors solver               " << endl;
	cout << "inverse a square matrix                       " << endl;
	cout << "all matrix are in cloumn major                " << endl;
	cout << "----------------------------------------------" << endl;
	bool cdtest = true;
	bool svdtest = true;
	bool eigentest = true;
	bool linearEqSolver = true;
	bool inverse = true;

	// cholesky decomposition
	if (cdtest) {

		cout << endl << endl;
		cout << "----------------------------------------------" << endl;
		cout << "Cholesky decomposition                        " << endl;
		cout << "----------------------------------------------" << endl;
		// basic setting for testing matrix
		UInt rowa = 3;
		UInt cola = 3;

		// set up matrix data
		// we need to set up a full rank matrix
		vector<Double> a(rowa*cola);
		a[0+0*rowa] =   4;
		a[1+0*rowa] =  12;
		a[2+0*rowa] = -16;
		a[1+1*rowa] =  37;
		a[2+1*rowa] = -43;
		a[2+2*rowa] =  98;
		cout << "matrix A, only lower part is provided" << endl;
		for(UInt i=0; i<rowa; i++) {
			for(UInt j=0; j<cola; j++) {
				printf("%-12.6f  ",a[i+j*rowa]); 
			}
			printf("\n");
		}

		// now in comparison
		vector<Double> b(a);
		mllt(&b.front(),rowa,rowa);
		vector<Double> c(rowa*cola);
		mmul(&b.front(),&b.front(),&c.front(),rowa,rowa,rowa,rowa,'N','T',ONE,ZERO);
		cout << "result matrix should be same with A" << endl;
		for(UInt i=0; i<rowa; i++) {
			for(UInt j=0; j<cola; j++) {
				printf("%-12.6f  ",c[i+j*rowa]); 
			}
			printf("\n");
		}
	}

	//testing SVD
	if (svdtest) {

		cout << endl << endl;
		cout << "----------------------------------------------" << endl;
		cout << "SVD                                           " << endl;
		cout << "----------------------------------------------" << endl;
		// define the data testing matrix
		UInt row = 3;
		UInt col = 2;
		vector<Double> x(row*col);
		UInt n = 1;
		for(UInt j=0; j<col; j++) {
			for(UInt i=0; i<row; i++) {
				x[i+j*row] = n;
				n++;
			}
		}
		for(UInt i=0; i<row; i++) {
			for(UInt j=0; j<col; j++) {
				printf("%-12.6f  ",x[i+j*row]); 
			}
			printf("\n");
		}

		vector<Double> S(row*row);
		vector<Double> VT(col*col);
		vector<Double> eigen(2);
		svd(&x.front(),&S.front(),&VT.front(),&eigen.front(),row,col,row);
		cout << "S and VT are all given from row 1 to row 2 to row 3 ..." << endl;
		cout << "S = -0.4287 0.806 0.4082 -0.5663 0.1124 -0.8165 -0.7039 -0.5812 0.4082 " 
			<< endl;
		cout << "result S:" << endl;
		for(UInt i=0; i<row; i++) {
			for(UInt j=0; j<row; j++) {
				printf("%-12.6f  ",S[i+j*row]); 
			}
			printf("\n");
		}
		cout << "VT = -0.3863 -0.9224 -0.9224 0.3863 " << endl;
		cout << "result V^{T}:" << endl;
		for(UInt i=0; i<col; i++) {
			for(UInt j=0; j<col; j++) {
				printf("%-12.6f  ",VT[i+j*col]); 
			}
			printf("\n");
		}
		cout << "eigen = 9.508 0.7729  " << endl; 
		cout << "result eigen value " << eigen[0] << " " << eigen[1] << endl;
	}

	// eigen test
	if (eigentest) {

		cout << endl << endl;
		cout << "----------------------------------------------" << endl;
		cout << "EIGEN decomposition                           " << endl;
		cout << "----------------------------------------------" << endl;
		// set up the original matrix
		UInt n = 5;
		vector<Double> b(n*n,ZERO);
		b[0+n*0] = 0.67;
		b[1+n*0] = -0.20;
		b[2+n*0] = 0.19;
		b[3+n*0] = -1.06;
		b[4+n*0] = 0.46;
		b[1+n*1] = 3.82;
		b[2+n*1] = -0.13;
		b[3+n*1] = 1.06;
		b[4+n*1] = -0.48;
		b[2+n*2] = 3.27;
		b[3+n*2] = 0.11;
		b[4+n*2] = 1.10;
		b[3+n*3] = 5.86;
		b[4+n*3] = -0.98;
		b[4+n*4] = 3.54;
		lowToUpper(n,b);
		cout << "original matrix for eigen calculation:" << endl;
		for(UInt i=0; i<n; i++) {
			for(UInt j=0; j<n; j++) {
				printf("%-12.6f  ",b[i+j*n]); 
			}
			printf("\n");
		}

		// real work
		vector<Double> e(n);
		cout << "Selected eigenvalues: 0.43   2.14   3.37 " << endl;
		cout << "Selected eigenvectors (stored columnwise) " << endl;
		cout << " -0.98  -0.01  -0.08 " << endl;  
		cout << "  0.01   0.02  -0.93 " << endl; 
		cout << "  0.04  -0.69  -0.07 " << endl; 
		cout << " -0.18   0.19   0.31 " << endl; 
		cout << "  0.07   0.69  -0.13 " << endl; 
		meigen(&b.front(),&e.front(),n);
		cout << "calculated eigenvalues:" << endl;
		for(UInt i=0; i<n; i++) cout << e[i] << endl;
		cout << "result eigenvectors:" << endl;
		for(UInt i=0; i<n; i++) {
			for(UInt j=0; j<n; j++) {
				printf("%-12.6f  ",b[i+j*n]); 
			}
			printf("\n");
		}
	}

	// test linear equation solver
	if(linearEqSolver) {

		cout << endl << endl;
		cout << "----------------------------------------------" << endl;
		cout << "AX = B solver                                 " << endl;
		cout << "----------------------------------------------" << endl;
		// matrix A
		// -5.86   3.99  -5.93  -2.82   7.69
		//  3.99   4.46   2.58   4.42   4.61
		// -5.93   2.58  -8.52   8.57   7.69
		// -2.82   4.42   8.57   3.72   8.07
		//  7.69   4.61   7.69   8.07   9.83
		UInt n = 5;
		vector<Double> A(n*n);
		A[0+n*0] = -5.86;
		A[1+n*0] =  3.99;
		A[2+n*0] = -5.93;
		A[3+n*0] = -2.82;
		A[4+n*0] =  7.69;
		A[1+n*1] =  4.46;
		A[2+n*1] =  2.58;
		A[3+n*1] =  4.42;
		A[4+n*1] =  4.61;
		A[2+n*2] = -8.52;
		A[3+n*2] =  8.57;
		A[4+n*2] =  7.69;
		A[3+n*3] =  3.72;
		A[4+n*3] =  8.07;
		A[4+n*4] =  9.83;

		// matrix B
		//    1.32  -6.33  -8.77
		//    2.22   1.69  -8.33
		//    0.12  -1.56   9.54
		//   -6.41  -9.49   9.56
		//    6.33  -3.67   7.48
		UInt rowB = n;
		UInt colB = 3;
		vector<Double> B(rowB*colB);
		B[0+n*0] =   1.32;
		B[1+n*0] =   2.22;
		B[2+n*0] =   0.12;
		B[3+n*0] =  -6.41;
		B[4+n*0] =   6.33;
		B[0+n*1] =  -6.33;
		B[1+n*1] =   1.69;
		B[2+n*1] =  -1.56;
		B[3+n*1] =  -9.49;
		B[4+n*1] =  -3.67;
		B[0+n*2] =  -8.77;
		B[1+n*2] =  -8.33;
		B[2+n*2] =   9.54;
		B[3+n*2] =   9.56;
		B[4+n*2] =   7.48;

		// DSYSV Example Program Results
		//
		// Solution
		//   1.17   0.52  -0.86
		//  -0.71   1.05  -4.90
		//  -0.63  -0.52   0.99
		//  -0.33   0.43   1.22
		//   0.83  -1.22   1.96
		//
		// Details of factorization
		//  -5.86   0.00   0.00   0.00   0.00
		//  -0.68   7.18   0.00   0.00   0.00
		//   1.01  -0.20  -2.82   0.00   0.00
		//   0.48   0.35  11.93   4.21   0.00
		//  -1.31   1.37   0.02   0.16   6.22
		printf("\n");
		printf("%s\n", "AX=B for symmetrical A");
		symAXeqB('L', &A.front(), rowB, &B.front(), colB);
		cout << " result B should be:" << endl;
		cout << "   1.17   0.52  -0.86 " << endl;
		cout << "  -0.71   1.05  -4.90 " << endl;
		cout << "  -0.63  -0.52   0.99 " << endl;
		cout << "  -0.33   0.43   1.22 " << endl;
		cout << "   0.83  -1.22   1.96 " << endl;
		cout << " result B after computation:" << endl;
		for(UInt i=0; i<rowB; i++) {
			for(UInt j=0; j<colB; j++) {
				printf("%-12.6f  ",B[i+j*rowB]); 
			}
			printf("\n");
		}
		cout << " result matrix A should be: " << endl;
		cout << "  -5.86   0.00   0.00   0.00   0.00 " << endl;
		cout << "  -0.68   7.18   0.00   0.00   0.00 " << endl;
		cout << "   1.01  -0.20  -2.82   0.00   0.00 " << endl;
		cout << "   0.48   0.35  11.93   4.21   0.00 " << endl;
		cout << "  -1.31   1.37   0.02   0.16   6.22 " << endl;
		cout << " result A after computation:" << endl;
		for(UInt i=0; i<n; i++) {
			for(UInt j=0; j<n; j++) {
				printf("%-12.6f  ",A[i+j*n]); 
			}
			printf("\n");
		}

		// DGESV Example Program Results
		// actually we use the same testing data
		// the factorization data is different from above
		// since it uses the LU decomposition
		//
		// Solution
		//   1.17   0.52  -0.86
		//  -0.71   1.05  -4.90
		//  -0.63  -0.52   0.99
		//  -0.33   0.43   1.22
		//   0.83  -1.22   1.96
		//
		A[0+n*0] = -5.86;
		A[1+n*0] =  3.99;
		A[2+n*0] = -5.93;
		A[3+n*0] = -2.82;
		A[4+n*0] =  7.69;
		A[1+n*1] =  4.46;
		A[2+n*1] =  2.58;
		A[3+n*1] =  4.42;
		A[4+n*1] =  4.61;
		A[2+n*2] = -8.52;
		A[3+n*2] =  8.57;
		A[4+n*2] =  7.69;
		A[3+n*3] =  3.72;
		A[4+n*3] =  8.07;
		A[4+n*4] =  9.83;
		B[0+n*0] =  1.32;
		B[1+n*0] =  2.22;
		B[2+n*0] =  0.12;
		B[3+n*0] = -6.41;
		B[4+n*0] =  6.33;
		B[0+n*1] = -6.33;
		B[1+n*1] =  1.69;
		B[2+n*1] = -1.56;
		B[3+n*1] = -9.49;
		B[4+n*1] = -3.67;
		B[0+n*2] = -8.77;
		B[1+n*2] = -8.33;
		B[2+n*2] =  9.54;
		B[3+n*2] =  9.56;
		B[4+n*2] =  7.48;
		lowToUpper(n,A);
		AXeqB(&A.front(), rowB, &B.front(), colB);
		printf("\n");
		printf("%s\n", "AX=B for asymmetrical A");
		cout << " result B should be:" << endl;
		cout << "   1.17   0.52  -0.86 " << endl;
		cout << "  -0.71   1.05  -4.90 " << endl;
		cout << "  -0.63  -0.52   0.99 " << endl;
		cout << "  -0.33   0.43   1.22 " << endl;
		cout << "   0.83  -1.22   1.96 " << endl;
		cout << " result B after computation:" << endl;
		for(UInt i=0; i<rowB; i++) {
			for(UInt j=0; j<colB; j++) {
				printf("%-12.6f  ",B[i+j*rowB]); 
			}
			printf("\n");
		}
	}

	if (inverse) {

		cout << endl << endl;
		cout << "----------------------------------------------" << endl;
		cout << "INVERSE MATRIX with LU decomposition          " << endl;
		cout << "----------------------------------------------" << endl;
		// input matrix
		UInt n = 3;
		vector<Double> A(n*n);
		A[0+0*n] =   4;
		A[1+0*n] =  12;
		A[2+0*n] = -16;
		A[1+1*n] =  37;
		A[2+1*n] = -43;
		A[2+2*n] =  98;
		lowToUpper(n,A);
		cout << "original matrix:" << endl;
		for(UInt i=0; i<n; i++) {
			for(UInt j=0; j<n; j++) {
				printf("%-12.6f  ",A[i+j*n]); 
			}
			printf("\n");
		}

		// now do inverse
		vector<Double> C(A);
		inverseGeneralMatrix(&A.front(),n);
		cout << "inverse of matrix:" << endl;
		for(UInt i=0; i<n; i++) {
			for(UInt j=0; j<n; j++) {
				printf("%-12.6f  ",A[i+j*n]); 
			}
			printf("\n");
		}

		// now test
		vector<Double> I(n*n);
		mmul(&A.front(),&C.front(),&I.front(),n,n,n,ONE,ZERO);
		cout << "should be identity matrix, with LU decomposition:" << endl;
		for(UInt i=0; i<n; i++) {
			for(UInt j=0; j<n; j++) {
				printf("%-12.6f  ",I[i+j*n]); 
			}
			printf("\n");
		}
	}

	if (inverse) {

		cout << endl << endl;
		cout << "----------------------------------------------" << endl;
		cout << "INVERSE MATRIX with CD decomposition          " << endl;
		cout << "----------------------------------------------" << endl;
		// input matrix
		UInt n = 4;
		vector<Double> A(n*n,ZERO);
		A[0+0*n] =  4.16;
		A[1+0*n] = -3.12;
		A[2+0*n] =  0.56;
		A[3+0*n] = -0.10;
		A[1+1*n] =  5.03;
		A[2+1*n] = -0.83;
		A[3+1*n] =  1.18;
		A[2+2*n] =  0.76;
		A[3+2*n] =  0.34; 
		A[3+3*n] =  1.18; 
		lowToUpper(n,A);
		cout << "original matrix, only lower triangular part is presented:" << endl;
		for(UInt i=0; i<n; i++) {
			for(UInt j=0; j<n; j++) {
				printf("%-12.6f  ",A[i+j*n]); 
			}
			printf("\n");
		}
		vector<Double> B(A);

		// now do inverse
		inverseSPDMatrix(&A.front(),n);
		lowToUpper(n,A);
		cout << "inverse of matrix:" << endl;
		for(UInt i=0; i<n; i++) {
			for(UInt j=0; j<n; j++) {
				printf("%-12.6f  ",A[i+j*n]); 
			}
			printf("\n");
		}

		// now test
		vector<Double> I(n*n,ZERO);
		mmul(&A.front(),&B.front(),&I.front(),n,n,n,ONE,ZERO);
		cout << "should be identity matrix, with CD decomposition:" << endl;
		for(UInt i=0; i<n; i++) {
			for(UInt j=0; j<n; j++) {
				printf("%-12.6f  ",I[i+j*n]); 
			}
			printf("\n");
		}
	}
}
