/**
 * CPP files corresponding to the block matrix class
 * \author fenglai liu 
 */
#include<string>
#include<cstdio>
#include "blas.h"
#include "blas1.h"
#include "excep.h"
#include "blockmatrix.h"
using namespace std;
using namespace blas;
using namespace excep;
using namespace blockmatrix;

void BlockMtrx::getData(const Mtrx& A, bool inTranspose) 
{
	if(inTranspose) {

		// for this case, the mother matrix's column is this
		// matrix row, we have to do it element by element
		for(UInt i=0; i<row; i++) {
			UInt c = i+rowPos;
			for(UInt j=0; j<col; j++) {
				UInt r = j+colPos;
				(*this)(i,j) += A(r,c);
			}
		}
	}else { 
		// for the non-transpose case, things is simple
		// we just use the copy from matrix function
		updateFromMatrix(rowPos,colPos,0,0,row,col,ONE,A);
	}
}

void BlockMtrx::updateData(Mtrx& A, bool inTranspose) const
{
	if(inTranspose) {

		// for this case, the mother matrix's column is this
		// matrix row, we have to do it element by element
		for(UInt i=0; i<row; i++) {
			UInt c = i+rowPos;
			for(UInt j=0; j<col; j++) {
				UInt r = j+colPos;
				A(r,c) += val(i,j); 
			}
		}
	}else { 
		// for the non-transpose case, things is simple
		// we just use the copy to matrix function
		updateToMatrix(0,0,rowPos,colPos,row,col,ONE,A);
	}
}

void BlockMtrx::add(const BlockMtrx& A) 
{
	if (same(A)) {
		UInt len = row*col;
		vaxpy(A.getPtr(),getPtr(),ONE,len);
	}else if (sameInTranspose(A)) {
		for(UInt i=0; i<row; i++) {
			for(UInt j=0; j<col; j++) {
				(*this)(i,j) += A(j,i); 
			}
		}
	}else{
		printf("this rowPos: %-12d  and the row: %-12d\n", (Int)rowPos, (Int)row);
		printf("this colPos: %-12d  and the col: %-12d\n", (Int)colPos, (Int)col);
		printf("A'   rowPos: %-12d  and the row: %-12d\n", (Int)A.getRowPos(), (Int)A.getRow());
		printf("A'   colPos: %-12d  and the col: %-12d\n", (Int)A.getColPos(), (Int)A.getCol());
		string infor = "The pass in block matrix data is not same/sameInTranspose as this one"; 
		Excep excep("BlockMtrx","add",EXCEPTION_BLOCK_MATRIX_INVALID_MATRIX_PASS_IN,infor);
		handleExcep(excep);
	}
}

void BlockMtrx::blockPrint(const string& title) const
{
	// now print out matrix content
	printf("%-15s %-15d\n","rowPos      : ", (Int)rowPos);
	printf("%-15s %-15d\n","colPos      : ", (Int)colPos);
	print(title,6);
}

