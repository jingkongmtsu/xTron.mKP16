//
// note:
//
// April 2014:
// transpose function are added and tested
//
// July 2014:
// rewrite the testing code for transpose
// test the in place and simple transpose timings and accuracy for both 
// single and double variables 
//
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<string> 
#include<vector>
#include "matrix.h"
#include "tbb/tbb.h"
using namespace matrix;
using namespace std;
using namespace tbb;

void para_matrix_test()
{
	// lapcck testing
	cout << "----------------------------------------------" << endl;
	cout << "normal matrix testing for parallel functions  " << endl;
	cout << "----------------------------------------------" << endl;
	cout << endl;
	bool testSelfOper          = true;

	// transpose function etc.
	if (testSelfOper) {

		cout << endl;
		cout << endl;
		cout << "----------------------------------------------" << endl;
		cout << "operation with itself                         " << endl;
		cout << "transpose operation test with timing          " << endl;
		cout << "here we actually do not test parallel         " << endl;
		cout << "what we really test here is transpose timing  " << endl;
		cout << "----------------------------------------------" << endl;

		for(UInt t=0; t<2; t++) {

			// set the in place transformation
			bool inPlace = false;
			if (t == 1) {
				cout << endl << endl;
				cout << "***********************" << endl;
				cout << "*  SIMPLE TRANSPOSE   *" << endl;
				cout << "***********************" << endl;
			}else{
				cout << "***********************" << endl;
				cout << "*  INPLACE TRANSPOSE  *" << endl;
				cout << "***********************" << endl;
				inPlace = true;
			}

			bool testSquare = true;
			if (testSquare) {

				// set up matrix
				UInt rowb = 30000;
				UInt colb = 30000;
				Mtrx b(rowb,colb);
				for(UInt i=0; i<colb; i++) {
					for(UInt j=0; j<rowb; j++) {
						b(j,i) = i+2*j+1;
					}
				}
				Mtrx a(b);

				// do transpose
				tick_count t0 = tick_count::now();
				b.transpose(inPlace);
				tick_count t1 = tick_count::now();
				Double t = (t1-t0).seconds();
				printf("%s  %-12.6f\n", "matrix transpose with with seconds for 30000*30000 matrix", t);

				// we test the accuracy
				Double thresh = 1.0E-5;
				b.transpose(false);
				b.diffComp(a,thresh,false);
			}

			bool testMNMatrix = true;
			if (testMNMatrix) {

				// set up matrix
				UInt rowb = 30000;
				UInt colb = 40000;
				Mtrx b(rowb,colb);
				for(UInt i=0; i<colb; i++) {
					for(UInt j=0; j<rowb; j++) {
						b(j,i) = i+2*j+1;
					}
				}
				Mtrx a(b);

				// do transpose
				tick_count t0 = tick_count::now();
				b.transpose(inPlace);
				tick_count t1 = tick_count::now();
				Double t = (t1-t0).seconds();
				printf("%s  %-12.6f\n", "matrix transpose with with seconds for 30000*40000 matrix", t);

				// we test the accuracy
				Double thresh = 1.0E-5;
				b.transpose(false);
				b.diffComp(a,thresh,false);
			}
		}
	}
}
