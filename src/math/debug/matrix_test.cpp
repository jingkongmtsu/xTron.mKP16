//
// note:
// April 2014, the following functions are well tested for normal matrix (single and double):
//
//	dimension change/memory handle test  
//	enlarge/shrink matrix test           
//	memory release test  
//
//	copy function test                  
//	copy full matrix test               
//	submatrix copy test                 
//	vector-matrix form and load test    
//
//	operation with itself                        
//	transpose operation test                     
//	copy lower triangular part to upper test   
//
//	BLAS in matrix                              
//	matrix add: A = 2*B+A                       
//	matrix scale: A = a*A                       
//	matrix set: A(i,j) = x                      
//	matrix multiplication: C=A*B+C, A is symm   
//	matrix multiplication: C=A*C+C              
//	matrix multiplication: C=A^T*B^T            
//	matrix multiplication: C=A*B, B is submatrix
//	matrix vector multiplication                
//	direct sum forming                          
//
//	LAPACK in matrix                           
//	matrix inverse with LU                     
//	matrix inverse with CD                     
//	matrix power function                      
//	matrix rank calculation with SVD           
//	matrix eigen is tested in power function   
//
//	April 2014:
//	add test for getRowColIndexFromPackIndex() function
//
//	June 2014:
//	re-test all of things after replace vector<Double> with 
//	DoubleVec (scalable allocator) from Intel
//
//	July 2014:
//	add in test for addTranspose
//	add in test for in place transpose and simple transpose
//	test the transpose cases for both single and double cases
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
#include "matrix.h"
using namespace matrix;
using namespace std;
using namespace omp_setting;

void matrix_test()
{
	// firstly, let's concentrate on serial mode
	// turn openmp off if the BLAS lib use openmp
	omp_turnoff();

	// lapcck testing
	cout << "----------------------------------------------" << endl;
	cout << "normal matrix testing                         " << endl;
	cout << "----------------------------------------------" << endl;
	cout << endl;
	bool dimensionChange       = true;
	bool testCopyFunction      = true;
	bool testSelfOper          = true;
	bool testLapack            = true;
	bool testBlas              = true;
	bool testUtil              = true;

	if (dimensionChange) {

		cout << endl;
		cout << endl;
		cout << "----------------------------------------------" << endl;
		cout << "dimension change/memory handle test           " << endl;
		cout << "enlarge/shrink matrix test                    " << endl;
		cout << "memory release test                           " << endl;
		cout << "----------------------------------------------" << endl;
		// enlarge and shrink 
		Mtrx a(10,10);
		cout << "Matrix A original size " << a.getRow() << " " << a.getCol() << endl;
		a.enlarge(3,3);
		cout << "after enlarge, Matrix A current size " << a.getRow() << 
			" " << a.getCol() << endl;
		a.reset(4,7,false);
		cout << "after reset, Matrix A current size " << a.getRow() << 
			" " << a.getCol() << endl;

		// test mem release
		Mtrx x(1000,1000);
		const DoubleVec& vec = x.getVec();
		cout << "Matrix X current size " << vec.size() << endl;
		x.memRelease();
		cout << "After memory release matrix X current size " << vec.size() << endl;
	}

	// copy functions
	if(testCopyFunction) {

		bool testFullCopy = true;
		bool testPartCopy = true;
		bool testVecCopy  = true;
		cout << endl;
		cout << endl;
		cout << "----------------------------------------------" << endl;
		cout << "copy function test                            " << endl;
		cout << "copy full matrix test                         " << endl;
		cout << "submatrix copy test                           " << endl;
		cout << "vector-matrix form and load test              " << endl;
		cout << "----------------------------------------------" << endl;

		if (testFullCopy) {

			// basic setting for testing matrix
			UInt rowa = 3;
			UInt cola = 3;

			// set up matrix data
			// we need to set up a full rank matrix
			Mtrx a(rowa,cola);
			for(UInt i=0; i<cola; i++) {
				for(UInt j=0; j<rowa; j++) {
					a(j,i) = 2*i+3*j+1;
				}
			}
			a.print("original data for full storage matrix a");

			// full copy between matrix and matrix
			Mtrx a1(rowa,cola);
			a1.copyMatrix(a);
			a1.print("copy full to full, should be same with matrix a");
		}

		// testing the copy between sub matrix
		// set up a bigger matrix
		if (testPartCopy) {

			UInt rowb = 6;
			UInt colb = 6;
			Mtrx b(rowb,colb);
			for(UInt i=0; i<colb; i++) {
				for(UInt j=0; j<rowb; j++) {
					b(j,i) = i+2*j+1;
				}
			}
			b.print("Original full matrix testing in copyToMatrix and copyFromMatrix");

			// try to fetch a copy from a diagonal part
			Mtrx t1(4,4);
			b.copyToMatrix(1,1,2,2,2,2,t1);
			t1.print("copy B's (1,1) to testing matrix(copy to)'s (2,2), size is 2*2");

			// copy data from a off-diagonal part
			Mtrx t3(3,3);
			b.copyToMatrix(4,2,1,0,2,2,t3);
			t3.print("copy B's (4,2) to testing matrix(copy to)'s (1,0), size is 2*2");

			// copy data from a off-diagonal part
			Mtrx t4(3,3);
			t4.copyFromMatrix(4,2,1,0,2,2,b);
			t3.print("copy B's (4,2) to testing matrix(copy from)'s (1,0), size is 2*2");
		}

		// test the copy between vector and matrix
		if (testVecCopy) {

			UInt rowc = 6;
			UInt colc = 6;
			Mtrx c(rowc,colc);
			for(UInt i=0; i<colc; i++) {
				for(UInt j=0; j<rowc; j++) {
					c(j,i) = i+j+1;
				}
			}
			c.print("Original full matrix testing in vector<->matrix copy");

			// form vector
			UInt len = (rowc*(1+rowc))/2;
			DoubleVec v(len);
			c.formSymmVec(v);

			// now load in vector
			Mtrx d(rowc,colc);
			d.loadSymmVec(v);
			d.print("should be same with original matrix");
		}
	}

	// transpose function etc.
	if (testSelfOper) {

		cout << endl;
		cout << endl;
		cout << "----------------------------------------------" << endl;
		cout << "operation with itself                         " << endl;
		cout << "transpose operation test                      " << endl;
		cout << "copy lower triangular part to upper test      " << endl;
		cout << "----------------------------------------------" << endl;

		//
		// since we can make sure that the simple transpose is correct
		// therefore we debug the code with simple transpose way
		//
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

			// square matrix
			UInt rowb = 6;
			UInt colb = 6;
			Mtrx b(rowb,colb);
			for(UInt i=0; i<colb; i++) {
				for(UInt j=0; j<rowb; j++) {
					b(j,i) = i+2*j+1;
				}
			}
			Mtrx a(b);
			b.transpose(inPlace);

			// we test the accuracy
			Double thresh = 1.0E-5;
			b.transpose(false);
			b.diffComp(a,thresh,false);

			rowb = 6;
			colb = 8;
			Mtrx c(rowb,colb);
			for(UInt i=0; i<colb; i++) {
				for(UInt j=0; j<rowb; j++) {
					c(j,i) = i+2*j+1;
				}
			}
			Mtrx d(c);
			c.transpose(inPlace);

			// we test the accuracy
			c.transpose(false);
			c.diffComp(d,thresh,false);
		}

		// copy low to upper
		UInt rowb = 6;
		UInt colb = 6;
		Mtrx b(rowb,colb);
		for(UInt i=0; i<colb; i++) {
			for(UInt j=0; j<rowb; j++) {
				b(j,i) = i+2*j+1;
			}
		}
		b.print("original matrix");
		b.copyLowToUpper();
		b.print("matrix after copy low to upper");

		// test the add in transpose
		for(UInt i=0; i<colb; i++) {
			for(UInt j=0; j<rowb; j++) {
				b(j,i) = i+2*j+1;
			}
		}
		b.addTranspose();
		Mtrx c(rowb,colb);
		for(UInt i=0; i<colb; i++) {
			for(UInt j=0; j<rowb; j++) {
				c(j,i) = i+2*j+1+j+2*i+1;
			}
		}
		b.diffComp(c,1.0E-5,false);
	}

	// test the direct sum
	if(testBlas) {

		cout << endl;
		cout << endl;
		cout << "----------------------------------------------" << endl;
		cout << " BLAS in matrix                               " << endl;
		cout << " matrix add: A = 2*B+A                        " << endl;
		cout << " matrix scale: A = a*A                        " << endl;
		cout << " matrix set: A(i,j) = x                       " << endl;
		cout << " matrix multiplication: C=A*B+C, A is symm    " << endl;
		cout << " matrix multiplication: C=A*C+C               " << endl;
		cout << " matrix multiplication: C=A^T*B^T             " << endl;
		cout << " matrix multiplication: C=A*B, B is submatrix " << endl;
		cout << " matrix vector multiplication                 " << endl;
		cout << " direct sum forming                           " << endl;
		cout << "----------------------------------------------" << endl;
		bool testadd          = true;
		bool testscale        = true;
		bool testSymmMult     = true;
		bool testSideMult     = true;
		bool testNormMult     = true;
		bool testSubMtrxMult  = true;
		bool testDirectSum    = true;

		// add function
		if(testadd) {
			Mtrx a(3,3);
			a.set(ONE);
			Mtrx b(3,3);
			b.set(TWO);
			a.add(b,TWO);
			a.print("element should be all five, testing add function");
		}

		// scale function
		if(testscale) {
			Mtrx a(3,3);
			a.set(ONE);
			a.scale(TWO);
			a.print("element should be all two, testing scale function");
		}

		// symm mult
		if (testSymmMult) {
			Mtrx x(3,3);
			x.set(ONE);
			Mtrx y(3,4);
			for(UInt i=0; i<4; i++) {
				for(UInt j=0; j<3; j++) {
					y(j,i) = i+2*j+1;
				}
			}
			Mtrx z(3,4);
			z.set(ONE);
			z.symMatrixMult(x,y,'L',ONE,ONE);
			y.print("y");
			z.print("x is all one, original Z is all one, final Z=X*Y+Z in symMatrixMult");
		}

		// side mult
		if (testSideMult) {
			Mtrx x(3,3);
			x.set(ONE);
			Mtrx y(3,4);
			for(UInt i=0; i<4; i++) {
				for(UInt j=0; j<3; j++) {
					y(j,i) = i+2*j+1;
				}
			}
			Mtrx z(y);
			z.print("original y, x is all one");
			y.sideMult(x,'L',false,ONE,ONE);
			y.print("final Y=X*Y+Y in sideMult, X pretend to be full matrix");
			y.copyMatrix(z);
			y.sideMult(x,'L',true,ONE,ONE);
			y.print("final Y=X*Y+Y in sideMult, X as symmetrical matrix");
		}

		// normal mult
		if (testNormMult) {
			Mtrx x(2,3);
			x.set(ONE);
			Mtrx y(3,4);
			for(UInt i=0; i<4; i++) {
				for(UInt j=0; j<3; j++) {
					y(j,i) = i+2*j+1;
				}
			}
			y.print("original y, x is all one");
			Mtrx z(4,2);
			z.mult(y,x,true,true,ONE,ZERO);
			z.print("z=y^t*x^t");
		}

		// sub matrix multiplication
		if (testSubMtrxMult) {

			// set the original matrix
			Mtrx x(10,10);
			x.set(ONE);
			Mtrx y(6,6);
			for(UInt i=0; i<6; i++) {
				for(UInt j=0; j<6; j++) {
					y(j,i) = i+2*j+1;
				}
			}
			Mtrx z(3,2);

			Mtrx smallx(3,3);
			smallx.copyFromMatrix(2,2,0,0,3,3,x);
			smallx.print("submatrix of x, from 2,2 in dimension 3*3, all element should be one");
			Mtrx smally(3,2);
			smally.copyFromMatrix(2,4,0,0,3,2,y);
			smally.print("submatrix of y, from 2,4 in dimension 3*2");
			y.print("this is original matrix");
			z.mult(x,y,2,2,2,4,3,3,3,2,false,false,ONE,ZERO);
			z.print("matrix multiplication in submatrix");
		}

		// direct sum
		if (testDirectSum) {

			Mtrx x(10,10);

			// piece matrix
			Mtrx t1(3,3);
			t1.set(ONE);
			Mtrx t2(3,3);
			t2.set(TWO);
			Mtrx t3(4,4);
			t3.set(THREE);

			// add in components
			UInt rowPos = 0;
			UInt colPos = 0;
			for(UInt iA=0; iA<3; iA++) {
				if (iA == 0) {
					x.addPieceInDirectSum(t1,rowPos,colPos);
					rowPos += t1.getRow();
					colPos += t1.getCol();
				}else if (iA == 1) {
					x.addPieceInDirectSum(t2,rowPos,colPos);
					rowPos += t2.getRow();
					colPos += t2.getCol();
				}else{
					x.addPieceInDirectSum(t3,rowPos,colPos);
					rowPos += t3.getRow();
					colPos += t3.getCol();
				}
			}
			x.print("matrix in direct sum test");
		}
	}

	// test the matrix operations
	if(testLapack) {

		cout << endl;
		cout << endl;
		cout << "----------------------------------------------" << endl;
		cout << " LAPACK in matrix                             " << endl;
		cout << " matrix inverse with LU                       " << endl;
		cout << " matrix inverse with CD                       " << endl;
		cout << " matrix power function                        " << endl;
		cout << " matrix rank calculation with SVD             " << endl;
		cout << " matrix eigen is tested in power function     " << endl;
		cout << "----------------------------------------------" << endl;

		// set up the original matrix
		UInt n = 5;
		Mtrx b(n,n);
		b(0,0) = 0.67;
		b(1,0) = -0.20;
		b(2,0) = 0.19;
		b(3,0) = -1.06;
		b(4,0) = 0.46;
		b(1,1) = 3.82;
		b(2,1) = -0.13;
		b(3,1) = 1.06;
		b(4,1) = -0.48;
		b(2,2) = 3.27;
		b(3,2) = 0.11;
		b(4,2) = 1.10;
		b(3,3) = 5.86;
		b(4,3) = -0.98;
		b(4,4) = 3.54;

		Mtrx b1(b);
		b1.copyLowToUpper();
		b1.print("original matrix for doing inverse/power");

		// general inverse
		bool matrixTranspose = false;
		Mtrx t1(b1);
		t1.inverseGenMatrix();
		Mtrx t2(n,n);
		t2.mult(b1,t1,matrixTranspose,matrixTranspose,ONE,ZERO);
		t2.print("Should be identity matrix after general inverse");

		// symmetrical inverse
		t1.copyMatrix(b);
		t1.inverseSymmMatrix();
		t1.copyLowToUpper();
		t2.mult(b1,t1,matrixTranspose,matrixTranspose,ONE,ZERO);
		t2.print("Should be identity matrix after SPD inverse");

		// power
		t1.copyMatrix(b);
		t1.power(-0.5E0, false, 1.0E-7);
		t2.mult(t1,t1,matrixTranspose,matrixTranspose,ONE,ZERO);
		t1.mult(t2,b1,matrixTranspose,matrixTranspose,ONE,ZERO);
		t1.print("Should be identity matrix after A^-1/2");

		//rank calculation
		Mtrx X(4,3);
		X.set(ONE);
		cout << "testing rank, X is all one, dimension is 4*3" << endl;
		UInt r1 = X.rank(1.0E-6);
		cout << "final rank is " << r1 << endl;
		Mtrx Y(3,4);
		X.set(ONE);
		cout << "testing rank, Y is all one, dimension is 3*4" << endl;
		UInt r2 = Y.rank(1.0E-6);
		cout << "final rank is " << r2 << endl;
	}	

	// test the matrix utility
	if(testUtil) {

		cout << endl;
		cout << endl;
		cout << "----------------------------------------------" << endl;
		cout << " matrix utility functions                     " << endl;
		cout << " getRowColIndexFromPackIndex()                " << endl;
		cout << "----------------------------------------------" << endl;

		UInt n = 1000;
		for(UInt times=1; times<10; times++) {

			// get the row/col index
			// make sure row>col
			UInt rp = n - times*11 + 6;
			UInt cp = n - times*12;

			// get index
			UInt index = rp+cp*(n-1)-(cp*(cp-1))/2;

			// now release new rp cp
			UInt nrp, ncp;
			getRowColIndexFromPackIndex(n,index,nrp,ncp);
			if (nrp != rp || ncp != cp) {
				cout << "getRowColIndexFromPackIndex is failed" << endl;
			}
		}
	}


}
