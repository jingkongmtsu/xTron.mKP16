//
//  note:
//
//  June, 2014:
//  test getData, updateData and constructor, init function; single and double
//
//  July, 2014:
//  add test for transpose and add, single and double
//
//  Sep, 2014:
//  update the test in according to the revise design of block matrix
//  add in test for getDataWithLowerTriAng function
//
//  Nov. 2014: 
//  consider to disable the getDataWithLowerTriAng, logic is too complicated
//  and in fact it could be replaced by something else
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<string> 
#include<vector>
#include "omp_setting.h"
#include "scalablevec.h"
#include "blockmatrix.h"
using namespace blockmatrix;
using namespace std;
using namespace omp_setting;

void blockmatrix_test()
{
	// lapcck testing
	cout << "----------------------------------------------" << endl;
	cout << "block matrix testing                          " << endl;
	cout << "----------------------------------------------" << endl;
	cout << endl;
	bool testSize   = true;
	bool getUpdate  = true;
	bool addTest    = true;

	// test the block matrix size
	if(testSize) {
		cout << "double var's size " << sizeof(Double) << endl;
		BlockMtrx A(10,10);
		cout << "block matrix A(10*10, uninitialized pos)'s size " << sizeof(A) << endl;
		BlockMtrx B(1000,500,5,6);
		cout << "block matrix B(1000*500, rowPos=5; colPos=6)'s size " << sizeof(B) << endl;
	}

	if (getUpdate) {

		cout << endl;
		cout << endl;
		cout << "----------------------------------------------" << endl;
		cout << "getting/updating data test                    " << endl;
		cout << "constructor test                              " << endl;
		cout << "init function test                            " << endl;
		cout << "getData function test                         " << endl;
		cout << "updateData function test                      " << endl;
		cout << "----------------------------------------------" << endl;

		// construct mother matrix
		UInt rowa = 5;
		UInt cola = 6;
		Mtrx a(rowa,cola);
		for(UInt i=0; i<cola; i++) {
			for(UInt j=0; j<rowa; j++) {
				a(j,i) = 2*i+3*j+1;
			}
		}

		// first construct, then initilize
		a.print("mother matrix for testing");
		BlockMtrx x(10,10);
		x.init(3,3,1,1);
		x.getData(a,false);
		x.blockPrint("matrix block in dimension of 3*3 from A(1,1)");
		Mtrx b1(a);
		x.set(ONE);
		x.updateData(b1,false);
		b1.print("update A with 3*3 matrix from A(1,1), submatrix is all ONE");
		cout << endl;
		cout << endl;

		a.print("mother matrix for testing");
		x.clear();
		x.init(4,4,1,1);
		x.getData(a,true);
		x.blockPrint("matrix block in dimension of 4*4 from A(1,1), in transpose form");
		Mtrx b2(a);
		for(UInt i=0; i<4; i++) {
			for(UInt j=0; j<4; j++) {
				if (j<=i) {
					x(j,i) = 0;
				}else{
					x(j,i) = 1;
				}
			}
		}
		x.updateData(b2,true);
		x.blockPrint("submatrix will be added in");
		b2.print("update A with 4*4 matrix from A(1,1) in transpose form");
		cout << endl;
		cout << endl;

		a.print("mother matrix for testing");
		BlockMtrx y(2,2,2,1);
		y.getData(a,true);
		y.blockPrint("matrix block in dimension of 2*2 from A(1,2), in transpose form");
		Mtrx b3(a);
		for(UInt i=0; i<2; i++) {
			for(UInt j=0; j<2; j++) {
				if (j<=i) {
					y(j,i) = 0;
				}else{
					y(j,i) = 1;
				}
			}
		}
		y.updateData(b3,true);
		y.blockPrint("submatrix will be added in");
		b3.print("update A with 2*2 matrix from A(1,2) in transpose form");
		cout << endl;
		cout << endl;

		a.print("mother matrix for testing");
		y.clear();
		y.init(3,3,1,0);
		y.getData(a,true);
		y.blockPrint("matrix block in dimension of 3*3 from A(0,1), in transpose form");
		Mtrx b4(a);
		for(UInt i=0; i<3; i++) {
			for(UInt j=0; j<3; j++) {
				if (j<=i) {
					y(j,i) = 0;
				}else{
					y(j,i) = 1;
				}
			}
		}
		y.updateData(b4,true);
		y.blockPrint("submatrix will be added in");
		b4.print("update A with 3*3 matrix from A(0,1) in transpose form");
		cout << endl;
		cout << endl;

		a.print("mother matrix for testing");
		y.clear();
		y.init(3,3,0,1);
		y.getData(a,false);
		y.blockPrint("matrix block in dimension of 3*3 from A(0,1)");
		Mtrx b5(a);
		for(UInt i=0; i<3; i++) {
			for(UInt j=0; j<3; j++) {
				if (j<=i) {
					y(j,i) = 0;
				}else{
					y(j,i) = 1;
				}
			}
		}
		y.updateData(b5,false);
		y.blockPrint("submatrix will be added in");
		b5.print("update A with 3*3 matrix from A(0,1)");
		cout << endl;
		cout << endl;
	}

	if(addTest) {
		cout << endl;
		cout << endl;
		cout << "----------------------------------------------" << endl;
		cout << "adding data test                              " << endl;
		cout << "----------------------------------------------" << endl;

		// first construct, then initilize
		BlockMtrx x(10,10);
		x.init(3,3,1,0);
		for(UInt i=0; i<3; i++) {
			for(UInt j=0; j<3; j++) {
				x(j,i) = i+2*j+1;
			}
		}
		x.blockPrint("matrix x");

		// another block matrix
		BlockMtrx y(3,3,0,1);
		for(UInt i=0; i<3; i++) {
			for(UInt j=0; j<3; j++) {
				y(j,i) = i+2*j+1;
			}
		}
		x.add(y);
		cout << "matrix y is in same content of x, however; with row/col on 0 1" << endl;
		x.blockPrint("matrix x after add in y: x = x + x^t");
	}
}
