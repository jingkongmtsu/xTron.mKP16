//
// note:
// to the date of Feb 2014, for netlib library:
// BLAS1 is fully tested for single and double;
// BLAS2 is fully tested for single and double;
// BLAS3 is fully tested for single and double;
//
// April 2014 (MKL version is 11.1, multi-threading version):
// BLAS lib are well tested for intel lib for single and double
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<string> 
#include<vector>
#include<cmath>
#include "libgen.h"
#include "blas.h"
#include "blas1.h"
#include "omp_setting.h"
using namespace blas;
using namespace std;
using namespace omp_setting;

void blas_test(bool testingBLAS1, bool testingBLAS2, bool testingBLAS3)
{
	// firstly, let's concentrate on serial mode
	// turn openmp off if the BLAS lib use openmp
	omp_turnoff();

	// blas1 testing
	if (testingBLAS1) {

		// basic setting for testing in BLAS1
		UInt lena = 100342;
		vector<Double> a(lena);
		for(UInt i=0; i<lena; i++) {
			a[i] = i+1;
		}
		Double thresh = THRESHOLD_MATH;

		cout << "----------------------------------------------" << endl;
		cout << "BLAS1 TESTING                                 " << endl;
		cout << "functions testing here:                       " << endl;
		cout << "vset, maxSearch, rmsd                         " << endl;
		cout << "vdot,vaxpy,vcopy,vscal,vmul,ddot3,v3add       " << endl;
		cout << "vaxyaddz, vmul3add                            " << endl;
		cout << "----------------------------------------------" << endl;
		vector<Double> v1(lena,ZERO);

		// vset
		vset(&v1.front(),ONE,lena);
		for(UInt i=0; i<lena; i++) {
			Double r = ONE;
			if (fabs(r-v1[i])>thresh) {
				cout << "in vset " << i << " value is wrong" <<endl;
			}
		}

		// maximum value search, should be 10000
		Double max = maxSearch(&a.front(),&v1.front(),lena);
		cout << "maxSearch " << max << " standard result is " << lena-1 <<  endl;

		// rmsd calculation
		vector<Double> vv(lena,TWO);
		Double rmsdVal = rmsd(&vv.front(),&v1.front(),lena);
		cout << "rmsd " << rmsdVal << " standard result is " << ONE << endl;

		// dot and norm should be same
		Double d1 = vdot(&a.front(),&v1.front(),lena);
		Double result = (1+lena)*lena/2;
		if (fabs(d1-result)>thresh) {
			cout << "vdot " << d1 << " standard result is " << result << " does not match" << endl;
		}

		// vaxpy
		cout << "v1 = 2*a + v1, original v1 is all one" << endl;
		vaxpy(&a.front(),&v1.front(),TWO,lena);
		for(UInt i=0; i<lena; i++) {
			Double r = 2*(i+1) + 1;
			if (fabs(r-v1[i])>thresh) {
				cout << "in vaxpy " << i << " value is wrong" <<endl;
			}
		}

		// vcopy
		cout << "vcopy " << endl;
		vcopy(&a.front(),&v1.front(),lena);
		for(UInt i=0; i<lena; i++) {
			Double r = i+1;
			if (fabs(r-v1[i])>thresh) {
				cout << "in vcopy " << i << " value is wrong" <<endl;
			}
		}

		// vscal
		cout << "vscal, a is THREE " << endl;
		vscal(&v1.front(),THREE,lena);
		for(UInt i=0; i<lena; i++) {
			Double r = 3*(i+1);
			if (fabs(r-v1[i])>thresh) {
				cout << "in vscal " << i << " value is wrong" <<endl;
			}
		}

		// vdiv
		vdiv(&a.front(),&a.front(),&v1.front(),lena);
		cout << "vdiv, should be all one" << endl;
		for(UInt i=0; i<lena; i++) {
			if (fabs(1.0E0-v1[i])>thresh) {
				cout << "in vscal " << i << " value is wrong" <<endl;
			}
		}

		// ddot3
		cout << "ddot3 " << endl;
		v1.assign(lena,ONE);
		Double d3 = vdot3(&v1.front(),&v1.front(),&v1.front(),lena);
		cout << "ddot3 " << d3 << " standard result should be " << lena << endl;

		// v3add
		cout << "v3add " << endl;
		v1.assign(lena,0);
		v3add(&a.front(),&a.front(),&a.front(),&v1.front(),lena);
		for(UInt i=0; i<lena; i++) {
			Double r = 3*(i+1);
			if (fabs(r-v1[i])>thresh) {
				cout << "in v3add " << i << " value is wrong" <<endl;
			}
		}

		//vaxyaddz 
		cout << "vaxyaddz " << endl;
		vector<Double> v2(lena,ONE);
		vset(&v1.front(),TWO,lena);
		vaxyaddz(&a.front(),&v1.front(),THREE,&v2.front(),lena);
		for(UInt i=0; i<lena; i++) {
			Double r = 3*(i+1)*2+1;
			if (fabs(r-v2[i])>thresh) {
				cout << "r: " << r << " our " << v2[i] << endl;
				cout << "in vaxyaddz " << i << " value is wrong" <<endl;
			}
		}

		//vmul3add
		cout << "vmul3add " << endl;
		vset(&v1.front(),TWO,lena);
		vset(&v2.front(),ONE,lena);
		vmul3add(&a.front(),&v1.front(),&v1.front(),&v2.front(),lena);
		for(UInt i=0; i<lena; i++) {
			Double r = (i+1)*2*2+1.0;
			if (fabs(r-v2[i])>thresh) {
				cout << "r: " << r << " our " << v2[i] << endl;
				cout << "in vmul3add " << i << " value is wrong" <<endl;
			}
		}
	}

	// blas2 testing
	if (testingBLAS2) {

		cout << "----------------------------------------------" << endl;
		cout << "blas2 testing                                 " << endl;
		cout << "all matrix are in cloumn major                " << endl;
		cout << "----------------------------------------------" << endl;

		// basic setting for testing matrix
		UInt rowa = 3;
		UInt cola = 3;

		// set up matrix and vector
		// we need to set up a full rank matrix
		vector<Double> a(rowa*cola);
		for(UInt i=0; i<cola; i++) {
			for(UInt j=0; j<rowa; j++) {
				a[j+i*rowa] = j + 2*i + 2;
			}
		}
		cout << "matrix A" << endl;
		for(UInt i=0; i<rowa; i++) {
			for(UInt j=0; j<cola; j++) {
				printf("%-12.6f  ",a[i+j*rowa]); 
			}
			printf("\n");
		}

		// set up two vector
		vector<Double> x;
		vector<Double> y;

		// symmetrical testing
		x.assign(rowa, ONE);
		y.assign(rowa, ZERO);
		vmsymv(&a.front(), rowa, rowa, &x.front(), &y.front(), ONE, ONE, false);
		cout << "x is all one, y is all zero " << endl;
		cout << "The result vector after A*x+y, only A's lower part used " << endl;
		for(UInt i=0; i<rowa; i++) cout << y[i] << endl;

		// general testing
		x.assign(rowa, ONE);
		y.assign(rowa, ZERO);
		vmgemv(&a.front(), 'N', rowa, rowa, rowa, &x.front(), &y.front(), ONE, ONE, false);
		cout << "x is all one, y is all zero " << endl;
		cout << "The result vector after A*x+y" << endl;
		for(UInt i=0; i<rowa; i++) cout << y[i] << endl;
	}


	// blas3 testing
	if (testingBLAS3) {

		cout << "----------------------------------------------" << endl;
		cout << "blas3 testing                                 " << endl;
		cout << "all matrix are in cloumn major                " << endl;
		cout << "----------------------------------------------" << endl;

		// basic setting for testing matrix
		UInt rowa = 3;
		UInt cola = 4;

		// set up matrix data
		// we need to set up a full rank matrix
		vector<Double> A(rowa*cola);
		for(UInt i=0; i<cola; i++) {
			for(UInt j=0; j<rowa; j++) {
				A[j+i*rowa] = i+j+1;
			}
		}
		cout << "matrix A is assymetrical:" << endl;
		for(UInt i=0; i<rowa; i++) {
			for(UInt j=0; j<cola; j++) {
				printf("%-12.6f ", A[i+j*rowa]); 
			}
			printf("\n");
		}

		//
		// we note that we do not need to test the msymmul with 
		// LD version
		//

		// first, B is symmetrical and multiply on the left side
		UInt rowb = rowa;
		UInt colb = rowa;
		vector<Double> B(rowb*colb,ONE);
		UInt rowc = rowb;
		UInt colc = cola;
		vector<Double> C(rowc*colc);
		cout << "matrix B is in dimension of " << rowb << "*" << colb << endl;
		cout << "it's symmetrical and all one " << endl;
		cout << "result matrix C is in dimension of " << rowc << "*" << colc << endl;
		msymmul(&B.front(),&A.front(),&C.front(),rowb,cola,'L', ONE, ZERO);
		cout << "testing L side symmetrical matrix multiplication: C = B*A " << endl;
		cout << "matrix C" << endl;
		for(UInt i=0; i<rowc; i++) {
			for(UInt j=0; j<colc; j++) {
				printf("%-12.6f  ", C[i+j*rowc]); 
			}
			printf("\n");
		}

		// second, B is symmetrical and multiply on the right side
		rowb = cola;
		colb = cola;
		B.assign(rowb*colb,ONE);
		rowc = rowa;
		colc = colb;
		C.assign(rowc*colc,ZERO);
		cout << "matrix B is in dimension of " << rowb << "*" << colb << endl;
		cout << "it's all one " << endl;
		cout << "result matrix C is in dimension of " << rowc << "*" << colc << endl;
		msymmul(&B.front(),&A.front(),&C.front(),rowb,rowa,'R', ONE, ZERO);
		cout << "testing R side symmetrical matrix multiplication: C = A*B " << endl;
		cout << "matrix C" << endl;
		for(UInt i=0; i<rowc; i++) {
			for(UInt j=0; j<colc; j++) {
				printf("%-12.6f  ", C[i+j*rowc]); 
			}
			printf("\n");
		}

		//
		// now test the mmul, without LD version
		//
		rowb = cola;
		colb = cola+1;
	   B.assign(rowb*colb,ONE);
		rowc = rowa;
		colc = colb;
		C.assign(rowc*colc,ONE);
		mmul(&A.front(),&B.front(),&C.front(),rowa,colb,cola,ONE,ONE);
		cout << "testing matrix multiplication: C = A*B+C " << endl;
		cout << "matrix B is in dimension of " << rowb << "*" << colb  << " all one "<< endl;
		cout << "matrix C is all one at first" << endl;
		for(UInt i=0; i<rowc; i++) {
			for(UInt j=0; j<colc; j++) {
				printf("%-12.6f  ", C[i+j*rowc]); 
			}
			printf("\n");
		}

		//
		// now test the mmul, without LD version
		// but with transpose
		//
		rowb = cola+2;
		colb = rowa;
	   B.assign(rowb*colb,ONE);
		rowc = cola;
		colc = rowb;
		C.assign(rowc*colc,ONE);
		mmul(&A.front(),&B.front(),&C.front(),rowa,cola,rowb,colb,'T','T',ONE,ZERO); 
		cout << "testing matrix multiplication: C = A^T*B^T " << endl;
		cout << "matrix B is in dimension of " << rowb << "*" << colb  << " all one "<< endl;
		cout << "matrix C" << endl;
		for(UInt i=0; i<rowc; i++) {
			for(UInt j=0; j<colc; j++) {
				printf("%-12.6f  ", C[i+j*rowc]); 
			}
			printf("\n");
		}

		//
		// finally, for the mmul with LD
		//
		cout << "the testing for mmul with leading dimension will be done in the other place" 
			<< endl;
	}
}
