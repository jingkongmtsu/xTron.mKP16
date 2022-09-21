/**
 * CPP files corresponding to the matrix class
 * \author fenglai liu 
 */
#include <boost/lexical_cast.hpp>
#include<cmath>
#include<cstdio>
#include<iostream>
#include<iomanip> 
#include "blas.h"
#include "blas1.h"
#include "lapack.h"
#include "matrix.h"
#include "inplacetranspose.h"
using namespace blas;
using namespace lapack;
using namespace matrix;

////////////////////////////////////////
//!!!  utility functions
////////////////////////////////////////
void matrix::getRowColIndexFromPackIndex(const UInt& n, const UInt& index, 
		UInt& rowIndex, UInt& colIndex)
{

	//
	// sometimes we may need to compute row/column index
	// for pack form data - that is to say, only the lower
	// triangular part plus the diagonal part is presented.
	//
	// therefore, this function get the row/column index
	// from the total index, which ranges (0-(1+n)*n/2)
	//
	
	//
	// keep an eye on n, it should be >=1
	//
#ifdef DEBUG
	if (n<1) {
		string infor = "matrix dimension should not be zero or negative";
		Excep excep("matrix","getRowColIndexFromPackIndex",
				EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}

	UInt lastIndex = n*(n+1)/2 - 1;
	if (index>lastIndex) {
		string infor = "the pass in index value is too large so it's invalid";
		Excep excep("matrix","getRowColIndexFromPackIndex",
				EXCEPTION_MATRIX_DIMENSION_OVERFLOW,infor);
		handleExcep(excep);
	}
#endif

	//
	// set the result to NULL state
	//
	rowIndex = -1;
	colIndex = -1;

	// 
	// firstly, is it the first column?
	//
	UInt inc = n-1;
	if (index<=inc) {
		rowIndex = index;
		colIndex = 0;
		return;
	}

	// get column index first
	for(UInt i=1; i<n; i++) {
		inc += n-i;
		if (index<=inc) {
			colIndex = i;
			break;
		}
	}

	// now get the row index
	rowIndex = index + (colIndex*(colIndex-1))/2 - colIndex*(n-1);
}

////////////////////////////////////////
//!!!  constructor section
////////////////////////////////////////
Mtrx::Mtrx(const Mtrx& A) {
	row = A.getRow();      
	col = A.getCol();
	// here we can not copy the leading dimension information
	// since we are unable to access to the ld of A
	ld  = row; 
	if (A.vecLen() > 0) {
		array.assign(row*col,ZERO);
		copyMatrix(A);
	}
}

Mtrx& Mtrx::operator=(const Mtrx& A) {
	if (this == &A) {
		return *this;
	}else{
		row = A.getRow();      
		col = A.getCol();
		// here we can not copy the leading dimension information
		// since we are unable to access to the ld of A
		ld  = row; 
		if (A.vecLen() > 0) {
			array.assign(row*col,ZERO);
			copyMatrix(A);
		}
		return *this;
	}
}

void Mtrx::reset(const UInt& nRow, const UInt& nCol, bool reform) 
{
	// test that whether we need reform
	// this happens when the len is larger than
	// the vector capacity
	bool needReform = reform;
	UInt len = nRow*nCol;
	if (len>array.capacity()) needReform = true;

	// reforming the data if necessary
	// we may raise a warning message if reformed is not required
	// but actually we need to do it
	if (needReform) {
		if (!reform && needReform) {
			string info = "reform matrix is not demanded, is there something wrong?"; 
			Excep excep("Mtrx","reset",EXCEPTION_MATRIX_RESET_WARNING,info);
			handleExcep(excep);
		}
		array.resize(len);
	}

	// reset the data
	row = nRow;
	col = nCol;
	ld  = row; 
	if (reform) {
		if (array.size()>=len) {
			vset(&array.front(),ZERO,len);
		}else{
			array.assign(len,ZERO);
		}
	}
}

void Mtrx::memRelease()
{
	// set up an empty vector
	// and swap their contents
	DoubleVec t;
	t.swap(array);
	row = 0;
	col = 0;
	ld  = 0;
}

void Mtrx::enlarge(const UInt& r, const UInt& c)
{
	// set up a tmp matrix to store old data
	Mtrx tmp(*this);

	// reform the data
	UInt newRow = row + r;
	UInt newCol = col + c;
	bool reform = true;
	reset(newRow,newCol,reform);

	// copy the old matrix back to the new one
	UInt rPosNewMtrx = 0;
	UInt cPosNewMtrx = 0;
	UInt rPosOldMtrx = 0;
	UInt cPosOldMtrx = 0;
	UInt rowOldMtrx  = tmp.getRow();
	UInt colOldMtrx  = tmp.getCol();
	copyFromMatrix(rPosOldMtrx,cPosOldMtrx,rPosNewMtrx,cPosNewMtrx, 
			rowOldMtrx,colOldMtrx,tmp); 
}

////////////////////////////////////////
//!!!  copy functions section
////////////////////////////////////////
void Mtrx::copyMatrix(const Mtrx& A) 
{
	// pre-check
#ifdef DEBUG
	if (row != A.getRow() || col != A.getCol()) {
		string infor = "this matrix's row and column should equal to A's dimension";
		Excep excep("Mtrx","copyMatrix",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif

	// for matrix is size of 0
	// just return
	if (A.vecLen() == 0) return; 

	if (A.getLd() == A.getRow() && ld == row) {
		vcopy(A.getPtr(),&array.front(),row*col);
	}else{
		// the matrix A is in part, that is to say,
		// leading dimension is not the row so we do it each column
		for(UInt j=0; j<col; j++) {
			const Double* p = A.getPtr(0,j);
			Double* ptr = getPtr(0,j);
			vcopy(p,ptr,row);
		}
	}
}

void Mtrx::copyFromMatrix(const UInt& rpos1, const UInt& cpos1, 
		const UInt& rpos2, const UInt& cpos2, const UInt& r, 
		const UInt& c, const Mtrx& A) 
{
	UInt readMatrixRowIndex  = rpos1;
	UInt writeMatrixRowIndex = rpos2;
	for(UInt j=0; j<c; j++) {
		UInt readMatrixColIndex = cpos1 + j;
		UInt writeMatrixColIndex = cpos2 + j;
		const Double* p = A.getPtr(readMatrixRowIndex,readMatrixColIndex);
		Double* ptr = getPtr(writeMatrixRowIndex,writeMatrixColIndex);
		vcopy(p,ptr,r);
	}
}

void Mtrx::copyToMatrix(const UInt& rpos1, const UInt& cpos1, 
		const UInt& rpos2, const UInt& cpos2, const UInt& r, 
		const UInt& c, Mtrx& A) const
{
	UInt readMatrixRowIndex  = rpos1;
	UInt writeMatrixRowIndex = rpos2;
	for(UInt j=0; j<c; j++) {
		UInt readMatrixColIndex = cpos1 + j;
		UInt writeMatrixColIndex = cpos2 + j;
		const Double* p = getPtr(readMatrixRowIndex,readMatrixColIndex);
		Double* ptr = A.getPtr(writeMatrixRowIndex,writeMatrixColIndex);
		vcopy(p,ptr,r);
	}
}

void Mtrx::updateFromMatrix(const UInt& rpos1, const UInt& cpos1, 
		const UInt& rpos2, const UInt& cpos2, const UInt& r, 
		const UInt& c, const Double& a, const Mtrx& A) 
{
	UInt readMatrixRowIndex  = rpos1;
	UInt writeMatrixRowIndex = rpos2;
	for(UInt j=0; j<c; j++) {
		UInt readMatrixColIndex = cpos1 + j;
		UInt writeMatrixColIndex = cpos2 + j;
		const Double* p = A.getPtr(readMatrixRowIndex,readMatrixColIndex);
		Double* ptr = getPtr(writeMatrixRowIndex,writeMatrixColIndex);
		vaxpy(p,ptr,a,r);
	}
}

void Mtrx::updateToMatrix(const UInt& rpos1, const UInt& cpos1, 
		const UInt& rpos2, const UInt& cpos2, const UInt& r, 
		const UInt& c, const Double& a, Mtrx& A) const
{
	UInt readMatrixRowIndex  = rpos1;
	UInt writeMatrixRowIndex = rpos2;
	for(UInt j=0; j<c; j++) {
		UInt readMatrixColIndex = cpos1 + j;
		UInt writeMatrixColIndex = cpos2 + j;
		const Double* p = getPtr(readMatrixRowIndex,readMatrixColIndex);
		Double* ptr = A.getPtr(writeMatrixRowIndex,writeMatrixColIndex);
		vaxpy(p,ptr,a,r);
	}
}

void Mtrx::copyFromVector(const UInt& rpos, const UInt& cpos, 
		const UInt& r, const UInt& c, const Double* v)
{
	UInt writeMatrixRowIndex = rpos;
	for(UInt j=0; j<c; j++) {
		const Double* p = v + j*r;
		UInt writeMatrixColIndex = cpos + j;
		Double* ptr = getPtr(writeMatrixRowIndex,writeMatrixColIndex);
		vcopy(p,ptr,r);
	}
}

void Mtrx::copyToVector(const UInt& rpos, const UInt& cpos, 
		const UInt& r, const UInt& c, Double* v) const
{
	UInt readMatrixRowIndex  = rpos;
	for(UInt j=0; j<c; j++) {
		UInt readMatrixColIndex = cpos + j;
		const Double* p = getPtr(readMatrixRowIndex,readMatrixColIndex);
		Double* ptr = v + j*r;
		vcopy(p,ptr,r);
	}

}

void Mtrx::updateFromVector(const UInt& rpos, const UInt& cpos, 
		const UInt& r, const UInt& c, const Double& a, const Double* v)
{
	UInt writeMatrixRowIndex = rpos;
	for(UInt j=0; j<c; j++) {
		UInt writeMatrixColIndex = cpos + j;
		const Double* p = v + j*r;
		Double* ptr = getPtr(writeMatrixRowIndex,writeMatrixColIndex);
		vaxpy(p,ptr,a,r);
	}
}

void Mtrx::updateToVector(const UInt& rpos, const UInt& cpos, 
		const UInt& r, const UInt& c, const Double& a, Double* v) const
{
	UInt readMatrixRowIndex  = rpos;
	for(UInt j=0; j<c; j++) {
		UInt readMatrixColIndex = cpos + j;
		const Double* p = getPtr(readMatrixRowIndex,readMatrixColIndex);
		Double* ptr = v + j*r;
		vaxpy(p,ptr,a,r);
	}
}

void Mtrx::formSymmVec(DoubleVec& vec) const {

	// additional check for the matrix
#ifdef DEBUG
	if (row != col) {
		string infor = "matrix dimension does not match";
		Excep excep("Mtrx","formSymmVec",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}

	// we also check the vector's dimension
	UInt len = (1+col)*col/2;
	if (vec.size() != len) {
		string infor = "the input vector's size is not matching matrix's dimension";
		Excep excep("Mtrx","formSymmVec",EXCEPTION_ARRAY_DIMENSION_NOT_MATCH,infor);
		handleExcep(excep);
	}
#endif

	// now produce the symm vec
	UInt pos = 0;
	for(UInt iCol=0; iCol<col; iCol++) {
		UInt len = col-iCol;
		vcopy(getPtr(iCol,iCol),&vec[pos],len);
		pos += len;
	}
}

void Mtrx::loadSymmVec(const DoubleVec& vec) {

	// additional check for the matrix
#ifdef DEBUG
	if (row != col) {
		string infor = "matrix is required to be square";
		Excep excep("Mtrx","loadSymmVec",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}

	// we also check the vector's dimension
	UInt len = (1+col)*col/2;
	if (vec.size() != len) {
		string infor = "the input vector's size is not matching matrix's dimension";
		Excep excep("Mtrx","loadSymmVec",EXCEPTION_ARRAY_DIMENSION_NOT_MATCH,infor);
		handleExcep(excep);
	}
#endif

	// now produce the symm vec
	UInt pos = 0;
	for(UInt iCol=0; iCol<col; iCol++) {
		UInt len = col-iCol;
		vaxpy(&vec[pos],getPtr(iCol,iCol),ONE,len);
		pos += len;
	}
}

////////////////////////////////////////
//!!!  transpose functions section
////////////////////////////////////////
void Mtrx::copyLowToUpper() {
#ifdef DEBUG
	if(row != col || row != ld){
		string infor = "square matrix is needed when you copy from lower->upper";
		Excep excep("Mtrx","copyLowToUpper",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif
	UInt n = row;
	for(UInt j=0; j<n; j++) {
		for(UInt i=j; i<n; i++) {
			array[j+i*n]=array[i+j*n];
		}
	}
}

void Mtrx::addTranspose() {
#ifdef DEBUG
	if(row != col || row != ld){
		string infor = "square matrix is needed when you do A = A + A^(T)";
		Excep excep("Mtrx","addTranspose",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif
	UInt n = row;
	for(UInt j=0; j<n; j++) {
		for(UInt i=j; i<n; i++) {
			array[i+j*n] += array[j+i*n];
			array[j+i*n]  = array[i+j*n];
		}
	}
}

void Mtrx::transposeInPlace()
{
#ifdef DEBUG
	if(row != ld) {
		string infor = "ld should equal to row in generalMatrixTranspose";
		Excep excep("Mtrx","transpose",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif

	// for square matrix, we can apply very simple technique
	// and this is also very fast
	// even for large matrix
	if (row == col) {
		for(UInt j=0; j<col; j++) {
			for(UInt i=j; i<row; i++) {
				Double tmp   = val(i,j);
				(*this)(i,j) = val(j,i);
				(*this)(j,i) = tmp;
			}
		}
		return;
	}

	// for the non-square matrices, sometimes we do not need to
	// do transpose - the transpose is just its self
	if (row<2 || col<2) return;

	// we will call a fortan code for the M*N matrix case
	// and M!=N
	Int M  = static_cast<Int>(row);
	Int N  = static_cast<Int>(col);
	Int MN = M*N;
	Int status = 0;

	// now let's evaluate the scratch size
	// we allocate the scratch on the stack
	// for the size of scratch, please read the 
	// comment of the fortran file 
	//
	// we note, that we do not check the status returned
	// it will never happen according to our arrangement
	// see the fortran code comment for more details
	//
	Int len = (M+N)/2;
	if (len<=100) {
		Int scr[100];
		Int lscr = 100;
		in_place_transpose(getPtr(),&M,&N,&MN,scr,&lscr,&status);
	}else if (len > 100 && len<=1000) {
		Int scr[1000];
		Int lscr = 1000;
		in_place_transpose(getPtr(),&M,&N,&MN,scr,&lscr,&status);
	}else if (len >1000 && len<=5000) {
		Int scr[5000];
		Int lscr = 5000;
		in_place_transpose(getPtr(),&M,&N,&MN,scr,&lscr,&status);
	}else if (len > 5000 && len<=10000) {
		Int scr[10000];
		Int lscr = 10000;
		in_place_transpose(getPtr(),&M,&N,&MN,scr,&lscr,&status);
	}else if (len > 10000 && len<=20000) {
		Int scr[20000];
		Int lscr = 20000;
		in_place_transpose(getPtr(),&M,&N,&MN,scr,&lscr,&status);
	}else if (len > 20000 && len<=30000) {
		Int scr[30000];
		Int lscr = 30000;
		in_place_transpose(getPtr(),&M,&N,&MN,scr,&lscr,&status);
	}else if (len > 30000 && len<=40000) {
		Int scr[40000];
		Int lscr = 40000;
		in_place_transpose(getPtr(),&M,&N,&MN,scr,&lscr,&status);
	}else if (len > 40000 && len<=50000) {
		Int scr[50000];
		Int lscr = 50000;
		in_place_transpose(getPtr(),&M,&N,&MN,scr,&lscr,&status);
	}else{
		string infor = "matrix dimension is too large for in place transpose";
		Excep excep("Mtrx","transposeInPlace",EXCEPTION_MATRIX_TRANSPOSE_ERROR,infor);
		handleExcep(excep);
	}

	// finally change the row and col
	UInt tmp = col;
	col = row;
	row = tmp;
	ld  = row;
}

void Mtrx::simpleTranspose()
{
#ifdef DEBUG
	if(row != ld) {
		string infor = "ld should equal to row in generalMatrixTranspose";
		Excep excep("Mtrx","transpose",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif

	// for the M*N matrix transpose, 
	// here we will resort to additional memory for a transpose
	DoubleVec vTmp(row*col);
	for(UInt i=0; i<col; i++) {
		for(UInt j=0; j<row; j++) {
			vTmp[i+j*col]=array[j+i*row];
		}
	}
	vcopy(&vTmp.front(),&array.front(),row*col);

	// finally change the row and col
	UInt tmp = col;
	col = row;
	row = tmp;
	ld  = row;
}

////////////////////////////////////////
//!!!  simple matrix operations 
////////////////////////////////////////
void Mtrx::scale(Double a) {
	if (ld == row) {
		UInt len = row*col;
		// in this case, the data is store continuously
		// this should be most cases in practical use
		vscal(&array.front(),a,len);
	}else{
		// if A is not storing data continuously,
		// we have to do job according to each column
		for(UInt i=0; i<col; i++) {
			Double* target = getPtr(0,i);
			vscal(target,a,row);
		}	
	}
}

void Mtrx::add(const Mtrx& M, Double a) 
{
	//
	// we note, that only if ld == row and M.getLd() == M.getRow()
	// the code is worthy of parallelized
	//

#ifdef DEBUG
	if (row!=M.getRow() || col!=M.getCol()) {
		string infor = "row/column should equal to input matrix's dimension";
		Excep excep("Mtrx","add",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif
	if (ld == row && M.getLd() == M.getRow()) {
		// in this case, the data is store continuously
		// this should be most cases in practical use
		UInt len = row*col;
		vaxpy(M.getPtr(),&array.front(),a,len);
	}else{
		// if either M or A is not storing data continuously,
		// we have to do job according to each column
		for(UInt i=0; i<col; i++) {
			const Double* source = M.getPtr(0,i);
			Double* target = getPtr(0,i);
			vaxpy(source,target,a,row);
		}	
	}
}

void Mtrx::set(Double a) {
	UInt len = row*col;
	array.assign(len,a);
}

Double Mtrx::rmsd(const Mtrx& A, bool withLowerPart) const 
{
#ifdef DEBUG
	// the input A and this matrix must be in same dimension
	if (row != A.getRow() || col != A.getCol()) {
		string infor = "matrix is required to be in same dimension of input matrix";
		Excep excep("Mtrx","rmsd",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif

	if (withLowerPart) {

#ifdef DEBUG
		// in this case the matrix must be square
		if (row != col) {
			string infor = "matrix is required to be square";
			Excep excep("Mtrx","rmsd",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
			handleExcep(excep);
		}
#endif

		// now loop over all of elements in lower part
		Double r = ZERO;
		for(UInt iCol=0; iCol<col; iCol++) {
			for(UInt iRow=iCol; iRow<row; iRow++) {
				Double v1 = val(iRow,iCol);
				Double v2 = A.val(iRow,iCol);
				Double v  = (v1-v2)*(v1-v2);
				r += v;
			}
		}

		// return the result
		UInt len = (row+1)*row/2;
		r = sqrt(r/len);
		return r;
	}

	// if we do not only use the lower part, then just call blas function
	UInt len = row*col;
	return blas::rmsd(getPtr(),A.getPtr(),len);
}

////////////////////////////////////////
//!!!  sophisticated matrix operations 
////////////////////////////////////////
void Mtrx::symMatrixMult(const Mtrx& A, const Mtrx& B, const char& state,
		Double alpha, Double beta, bool withPara)
{
#ifdef DEBUG
	if(A.getRow() != A.getCol()) {
		string infor = "matrix should be symmetrical";
		Excep excep("Mtrx","symMatrixMult",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
	if (state != 'R' && state != 'L') {
		string infor = "input char parameter is invalid!"; 
		Excep excep("Mtrx","symMatrixMult",EXCEPTION_CHAR_VARIABLE_INVALID,infor);
		handleExcep(excep);
	}
#endif
	UInt N = A.getRow();
	UInt M = B.getCol();
	if (state == 'R') M = B.getRow();
	msymmul(A.getPtr(),B.getPtr(),getPtr(),N,M,A.getLd(),B.getLd(),getLd(), 
			state,alpha,beta,withPara);
}

void Mtrx::sideMult(const Mtrx& A, const char& state, bool isSymmA, 
		Double alpha, Double beta, bool withPara) 
{
#ifdef DEBUG
	// check the input state
	if (state != 'R' && state != 'L') {
		string infor = "input char parameter is invalid!"; 
		Excep excep("Mtrx","sideMult",EXCEPTION_CHAR_VARIABLE_INVALID,infor);
		handleExcep(excep);
	}
#endif
	if (isSymmA) {

		// set dimension
		UInt N = A.getRow();
		UInt M = getCol();
		if (state == 'R') M = getRow();

		// set up a tmp matrix to store a copy of C
		Mtrx tmp(*this);
		msymmul(A.getPtr(),tmp.getPtr(),getPtr(),N,M,A.getLd(),tmp.getLd(),getLd(), 
				state,alpha,beta,withPara);

	}else{

		// set up a tmp matrix to store a copy of C
		Mtrx tmp(*this);
		if (state == 'R') {
			// C = C*A + C
			mmul(tmp.getPtr(),A.getPtr(),getPtr(),tmp.getRow(),tmp.getCol(),A.getRow(),
					A.getCol(),tmp.getLd(),A.getLd(),ld,'N','N',alpha,beta,withPara);
		}else if (state == 'L') {
			// C = A*C + C
			mmul(A.getPtr(),tmp.getPtr(),getPtr(),A.getRow(),
					A.getCol(),tmp.getRow(),tmp.getCol(),A.getLd(),tmp.getLd(),ld,
					'N','N',alpha,beta,withPara);
		}
	}
}

void Mtrx::mult(const Mtrx& A, const Mtrx& B, bool transposeA, bool transposeB, 
		Double alpha, Double beta,bool withPara) 
{
	// check the matrix status
	char symA = 'N';
	if (transposeA) symA = 'T';
	char symB = 'N';
	if (transposeB) symB = 'T';
	mmul(A.getPtr(),B.getPtr(),&array.front(),A.getRow(),A.getCol(),B.getRow(),
			B.getCol(),A.getLd(),B.getLd(),ld,symA,symB,alpha,beta,withPara);
}

void Mtrx::mult(const Mtrx& A, const Mtrx& B, const UInt& rowAPos, 
		const UInt& colAPos, const UInt& rowBPos, const UInt& colBPos,
		const UInt& rowA, const UInt& colA, const UInt& rowB, const UInt& colB, 
		bool transposeA, bool transposeB, Double alpha, Double beta, bool withPara)
{
#ifdef DEBUG
	if (rowAPos >=A.getRow() || colAPos>=A.getCol() || 
			rowBPos>=B.getRow() || colBPos>=B.getCol()) {
		string infor  = "please check the position for row/column of A or B";
		Excep excep("Mtrx","mult for submatrix",EXCEPTION_MATRIX_DIMENSION_OVERFLOW,infor);
		handleExcep(excep);
	}
	if (rowA >A.getRow() || colA>A.getCol() || 
			rowB>B.getRow() || colB>B.getCol()) {
		string infor  = "please check the number of row/column for A or B";
		Excep excep("Mtrx","mult for submatrix",EXCEPTION_MATRIX_DIMENSION_OVERFLOW,infor);
		handleExcep(excep);
	}
#endif

	// check the matrix status
	char symA = 'N';
	if (transposeA) symA = 'T';
	char symB = 'N';
	if (transposeB) symB = 'T';

	// real work
	const Double* aptr = A.getPtr(rowAPos,colAPos);
	const Double* bptr = B.getPtr(rowBPos,colBPos);
	mmul(aptr,bptr,&array.front(),rowA,colA,rowB,colB,A.getLd(),
			B.getLd(),ld,symA,symB,alpha,beta,withPara);
}

void Mtrx::mmulv(const Double* x, Double* y, const Double& alpha, 
		const Double& beta, bool withLowerTriAngularPart, 
		bool needTranspose, bool withPara) const 
{
	if (withLowerTriAngularPart) {
		vmsymv(&array.front(),getRow(),getLd(),x,y,alpha,beta,withPara);
	}else{
		char symA = 'N';
		if (needTranspose) symA = 'T';
		vmgemv(&array.front(), symA, getRow(), getCol(),getLd(),x,y,alpha,beta,withPara);
	}
}

void Mtrx::addPieceInDirectSum(const Mtrx& A, const UInt& rowPos, const UInt& colPos) 
{
#ifdef DEBUG
	if (row != ld) {
		string infor = "matrix leading dimension should equal to row";
		Excep excep("Mtrx","addPieceInDirectSum",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
	if (row != col) {
		string infor = "matrix must be square";
		Excep excep("Mtrx","addPieceInDirectSum",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif
	if (rowPos + A.getRow() > row) {
		string infor1 = "row limit is: " + boost::lexical_cast<string>(row); 
		string infor2 = " currrnt row will be: " + boost::lexical_cast<string>(rowPos+A.getRow());
		string infor  = infor1 + infor2;
		Excep excep("Mtrx","addPieceInDirectSum",EXCEPTION_MATRIX_DIMENSION_OVERFLOW,infor);
		handleExcep(excep);
	}
	if (colPos + A.getCol() > col) {
		string infor1 = "col limit is: " + boost::lexical_cast<string>(col); 
		string infor2 = " current col will be: " + boost::lexical_cast<string>(colPos+A.getCol());
		string infor  = infor1 + infor2;
		Excep excep("Mtrx","addPieceInDirectSum",EXCEPTION_MATRIX_DIMENSION_OVERFLOW,infor);
		handleExcep(excep);
	}
	UInt nRows = A.getRow();
	UInt nCols = A.getCol();
	for(UInt i=0; i<nCols; i++) {
		Double* ptr = getPtr(rowPos,colPos+i);
		const Double* source = A.getPtr(0,i);
		vcopy(source,ptr,nRows);
	}
}

Double Mtrx::dotProduct(const Mtrx& D, bool isSymm) const
{
	if (isSymm) {

		// additional check on dimension 
#ifdef DEBUG
		if (row != D.getRow() && col != D.getCol()) {
			string infor = "two matrices dimension should math with each other";
			Excep excep("Mtrx","dotProduct",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
			handleExcep(excep);
		}
#endif

		// for symmetry data matrix, we only use its lower triangular part
		// firstly, let's go to see the lower triangular part
		Double eLowTriAngular = ZERO;
		for(UInt iCol=0; iCol<getCol(); iCol++) {
			const Double* dataA = getPtr(iCol,iCol);
			const Double* dataD = D.getPtr(iCol,iCol);
			UInt len = getCol() - iCol;
			eLowTriAngular += vdot(dataA,dataD,len);
		}

		// now it's the diagonal part
		Double eDiagonal = ZERO;
		for(UInt iCol=0; iCol<getCol(); iCol++) {
			eDiagonal += val(iCol,iCol)*D(iCol,iCol);
		}

		// finally, get the whole value
		return eLowTriAngular*TWO - eDiagonal;
	
	}else{

		// for assymetrical matrix case, just multiply all of elements
		Double e = ZERO;
		for(UInt iCol=0; iCol<col; iCol++) {
			const Double* dataA = getPtr(0,iCol);
			const Double* dataD = D.getPtr(0,iCol);
			e += vdot(dataA,dataD,col);
		}
		return e;
	}
}

void Mtrx::inverseSymmMatrix(bool withPara) {
#ifdef DEBUG
	if (row != col) {
		string infor = "the matrix must be square";
		Excep excep("Mtrx","inverse",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif
	inverseSPDMatrix(&array.front(),row, withPara); 
}

void Mtrx::inverseGenMatrix(bool withPara) {
#ifdef DEBUG
	if (row != col) {
		string infor = "the matrix must be square";
		Excep excep("Mtrx","inverse",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif
	inverseGeneralMatrix(&array.front(),row,withPara); 
}

void Mtrx::power(Double p, bool withPara, Double linearDepThresh) 
{
	if (row != col) {
		string infor = "the matrix must be square";
		Excep excep("Mtrx","power",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}

	// the general idea here is to decompose the matrix
	// into the form of A = PDP^T
	// D is the eigen values
	// P stores the corresponding eigen vector
	// we require that the given matrix is symmetrical 
	//
	// also for the linear dependent case, we will rise an error
	// for this case

	// prepare work
	UInt n = row;
	DoubleVec D(n, ZERO);

	// get the eigen vector and eigen values for this symmetrical matrix
	Mtrx P(n,n);
	P.copyMatrix(*this);
	P.doEigen(D,withPara);

	// checking that the conditions
	// all of eigen values should be >=0 within the error range
	Double eMin = LARGE_FLOAT_NUMBER;
	for(UInt i=0; i<n; i++) {
		if (D[i] < eMin) eMin = D[i];
	}
	if (eMin<ZERO) {
		string infor = "integrals matrix calculation could be inaccurate, "
			"please increase integral calculation threshold";
		Excep excep("Mtrx","power",EXCEPTION_EIGEN_VALUE_LESS_THAN_ZERO,infor);
		handleExcep(excep);
	}

	// let's see whether we have linear dependency in the original matrix
	if (eMin<linearDepThresh && p<ZERO) {
		cout << "linear dependent threshold is " << linearDepThresh << endl;
		cout << "the minimum eigenvalue is " << eMin << endl;
		string infor = "linear dependency exists so we can not do the power function";
		Excep excep("Mtrx","power",EXCEPTION_LINEAR_DEPENDENCY_ERROR_IN_MATRIX,infor);
		handleExcep(excep);
	}

	// keep a copy of P
	Mtrx P1(P);

	// scale the P
	for(UInt i=0; i<n; i++) {
		Double f = pow(D[i],p);
		vscal(P.getPtr(0,i), f, n);
	}

	// form final matrix
	mult(P,P1,false,true,ONE,ZERO,withPara);
}

void Mtrx::doEigen(DoubleVec& eigenValue, bool withPara) 
{
#ifdef DEBUG
	if (row != col) {
		string infor = "the matrix must be square to do eigen decomposition";
		Excep excep("Mtrx","doEigen",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif

	//real work
	meigen(&array.front(),&eigenValue.front(),row,withPara);
}

void Mtrx::formOrthoMatrix(const Mtrx& ov, Double linearDepThresh, bool withPara) 
{
#ifdef DEBUG
	if (row != col) {
		string infor = "the matrix must be square";
		Excep excep("Mtrx","formOrthoMatrix",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
#endif

	//
	// the general idea here is two ways:
	//
	// 1 if the given matrix has eigenvalue near zero, that is to say;
	// it's singular matrix; then we return the linear independent 
	// eigen vectors and scale it with eigenvalue^(-0.5). This is
	// canonical orthogonalization procedure 
	//
	// 2 if the matrix is not singular, then we return S^{-1/2},
	// which is the lowdin orthogonalization procedure
	// in general the lowdin orthogonalization procedure is better
	// however, it's only can be applied to the linear independent matrix
	//
	// one thing to note. That here we do not check rank first.
	// Since SVD decomposition is much more expensive than eigen
	// decomposition. Therefore we only use eigen decomposition
	// to check the rank, like what we did in power function
	//

	// prepare work
	UInt n = row;
	DoubleVec D(n, ZERO);

	// get the eigen vector and eigen values for oringinal matrix
	// we note that only lower triangular part is needed
	copyMatrix(ov);
	doEigen(D,withPara);
	//print("after diagonalize");

	// checking that the conditions
	// all of eigen values should be >=0 within the error range
	Double eMin = LARGE_FLOAT_NUMBER;
	for(UInt i=0; i<n; i++) {
		if (D[i] < eMin) eMin = D[i];
	}
	if (eMin<ZERO) {
		string infor = "integrals matrix calculation could be inaccurate";
		Excep excep("Mtrx","formOrthoMatrix",EXCEPTION_EIGEN_VALUE_LESS_THAN_ZERO,infor);
		handleExcep(excep);
	}

	// let's see whether we have linear dependency in the overlap matrix
	bool hasLinearDependency = false;
	if (eMin<linearDepThresh) {
		hasLinearDependency = true;
		string infor = "linear dependency exists in the molecule of the given basis set";
		Excep excep("Mtrx","formOrthoMatrix",EXCEPTION_LINEAR_DEPENDENCY_IN_MATRIX,infor);
		handleExcep(excep);
	}

	// now let's process the normal Lowdin matrix
	if (! hasLinearDependency) {

		// we need P to store the result of S^{-1/2} 
		// also P1 to store the original eigenvectors
		Mtrx P(n,n);
		Mtrx P1(*this);

		// scale the P
		for(UInt i=0; i<n; i++) {
			Double f = pow(D[i],-0.5E0);
			vscal(getPtr(0,i), f, n);
		}

		// form final matrix
		P.mult(*this,P1,false,true,ONE,ZERO,withPara);
		copyMatrix(P);

		// now let's exit
		return;
	}

	// we will count the eigen vectors actually used 
	// will eliminate the eigen velues which are close to zero
	// copy the eigen-vectors whose eigen-value are large enough
	// we also scale the eigen vector here
	Mtrx P1(n,n);
	UInt colIndex = 0;
	for(UInt i=0; i<n; i++) {
		if (D[i] < linearDepThresh) {
			printf("%dth eigenvector is removed, and it's eigenvalue: %-16.10f, "
					"linear dependency threshold is %-16.14f\n", (Int)i, D[i], linearDepThresh);
		}else{
			Double f = pow(D[i],-0.5E0);
			vscal(getPtr(0,i), f, n);
			vcopy(getPtr(0,i),P1.getPtr(0,colIndex),n);
			colIndex++;
		}
	}
	if (colIndex<n) {
		string infor = boost::lexical_cast<string>(n-colIndex) + " vectors are removed";
		Excep excep("Mtrx","formOrthoMatrix",EXCEPTION_LINEAR_DEPENDENCY_IN_MATRIX,infor);
		handleExcep(excep);
	}else{
		printf("the input linear dependence threshold is %f\n", linearDepThresh);
		string infor = "all of matrix eigenvalues are less than the linear dependence threshold?";
		Excep excep("Mtrx","formOrthoMatrix",EXCEPTION_LINEAR_DEPENDENCY_ERROR_IN_MATRIX,infor);
		handleExcep(excep);
	}	
	UInt nDep = colIndex;
	P1.reset(n,nDep,false);
	//P1.print("result ortho matrix in linear dependence situation");

	// finally reset the matrix
	reset(n,nDep,true);
	copyMatrix(P1);

	// we need to test that whether the canonical orthogonalization procedure
	// is able to give good result
	// S*ortho (S could be only lower part matrix)
	P1.symMatrixMult(ov,*this,'L',ONE,ZERO);
	//P1.print("S*ortho");

	// ortho^{T}*S*ortho
	// this should be identity matrix
	Mtrx tmp(nDep,nDep);
	tmp.mult(*this,P1,true,false,ONE,ZERO);
	//tmp.print("ortho^{T}*S*ortho");

	// now let's check
	for(UInt i=0; i<nDep; i++) tmp(i,i) += MINUS_ONE;
	Double dev = maxSearch(tmp.getPtr(),nDep*nDep);
	if (dev>linearDepThresh) {
		tmp.print("the matrix of ortho^{T}*S*ortho after canonical orthogonalization procedure");
		string infor = "canonical orthogonalization procedure does not give a good result";
		Excep excep("Mtrx","formOrthoMatrix",EXCEPTION_LINEAR_DEPENDENCY_ERROR_IN_MATRIX,infor);
		handleExcep(excep);
	}
}

UInt Mtrx::rank(Double thresh, bool withPara) const
{
	// the general idea here is to use SVD to decompose 
	// the matrix. the input matrix should be symmetrical
	// matrix

	// prepare work
	// since SVD may broke the content of original matrix
	// here we need a copy
	Mtrx X(*this);
	Mtrx U(row,row);
	Mtrx VT(col,col);
	UInt len=row;
	if (col<row) len = col;
	DoubleVec S(len,ZERO);

	// now do SVD
	svd(X.getPtr(),U.getPtr(),VT.getPtr(),&S.front(),row,col,ld,withPara); 

	// let's see whether the eigenvalue are all > 0?
	// if it's not this case, we need to report an error
	bool allPositive = true;
	for(UInt i=0; i<len; i++) {
		if (S[i] < ZERO && fabs(S[i])>thresh) allPositive = false;
	}
	if (! allPositive) {
		string infor = "the matrix's eigenvalue is significantly less than zero, please check the reason";
		Excep excep("Mtrx","rank",EXCEPTION_EIGEN_VALUE_LESS_THAN_ZERO,infor);
		handleExcep(excep);
	}

	// since the eigenvalues are arranged in ascending way, therefore
	// we will count the range of eigen vectors actually used - this 
	// determines that whether the matrix is singular or not
	UInt count = 0;
	for(UInt i=0; i<len; i++) {
		if (S[i] < thresh) {
			count++;
		}
	}

	// now return the infor
	if (count>0) {
		string infor = boost::lexical_cast<string>(count) + " vectors are linear dependent";
		Excep excep("Mtrx","rank",EXCEPTION_LINEAR_DEPENDENCY_IN_MATRIX,infor);
		handleExcep(excep);
	}
	return len - count;
}

////////////////////////////////////////
//!!!  additional functions
////////////////////////////////////////
void Mtrx::print(string title, UInt cols) const 
{

	// get the number of rows and columns
	UInt rows   = getRow();
	UInt column = getCol();

	// check the number of rows
	// we do not do it 
	UInt iwidth = 5; // this is related to the order of rows
	string space = "     ";
	if (rows*column > 1000000) {
		string infor = "too many data required to printed out, just return";
		Excep excep("Mtrx","print",EXCEPTION_TOO_MANY_TO_PRINT_MATRIX,infor);
		handleExcep(excep);
		return;
	}

	// choose the printing precision in accordance with the column lines
	UInt width = 0; 
	UInt ndigits = 0; 
	switch (cols) {
		case 3:
			width   = 20; 
			ndigits = 10;
			break;
		case 4:
			width   = 16; 
			ndigits = 10;
			break;
		case 5:
			width   = 14; 
			ndigits = 8;
			break;
			break;
		case 6:
			width   = 12; 
			ndigits = 7;
			break;
		case 7:
			width   = 10; 
			ndigits = 5;
			break;
		default:
			string infor = "the width is set to 10, n digits is set to 5";
			Excep excep("Mtrx","print",EXCEPTION_INVALID_COLS_PRINT_MATRIX,infor);
			handleExcep(excep);
			width   = 10; 
			ndigits = 5;
			break;
	}
	cout.precision(ndigits);

	// real printing work
	UInt start = 0;
	UInt colIndex = 1;
	UInt end = 0;
	cout << title << endl;
	while(start<column) {

		// how much data we will print each time
		if (start+cols <= column) {
			end = start+cols;
		}else{
			end = column;
		}

		// print out the line number for the data
		cout << space << " ";
		for(UInt i=start; i<end; i++) {
			cout << setw(width) << colIndex << " ";
			colIndex++;
		}
		cout << endl;

		// print out the data
		for(UInt i=0; i<rows; i++) {
			cout << setw(iwidth) << left << i+1 << " ";
			for(UInt j=start; j<end; j++) {
				Double v = val(i,j);
				cout<<setw(width)<<showpoint<<fixed<<right<< v << " ";
			}
			cout << endl;
		}
		start = end;
	}
}

bool Mtrx::diffComp(const Mtrx& A, const Double& thresh, bool onlyCompareLowPart, 
		bool printMore) const
{
	Mtrx C(*this);
	C.add(A,-1.0E0);
	bool fail = false;
	Double largestDiff = ZERO;
	for(UInt j=0; j<C.getCol(); j++) {
		for(UInt i=0; i<C.getRow(); i++) {
			if (onlyCompareLowPart && i<j) continue;
			Double diff = fabs(C(i,j));
			if (diff > thresh) {
				fail = true;
				if (diff>largestDiff) largestDiff = diff;
				if (printMore) {
					Double per        = 1000.0E0;
					Double absAij     = fabs(val(i,j));
					if (absAij>1.0E-14) {
						per = diff/absAij;
					}
					printf("Element (%d , %d) the abs difference is %-16.10f and abs_diff/(i,j) is %-16.10f\n",
							(Int)i, (Int)j, diff, per);
				}
			}
		}
	}
	if (fail) {
		if (onlyCompareLowPart) {
			C.copyLowToUpper();
		}
		const DoubleVec& v = C.getVec();
		Double rmsdErr = blas::rmsd(&v[0],C.getRow()*C.getCol());
		printf("maximum error is %-16.10f\n",largestDiff);
		printf("RMSD error is %-16.10f\n",rmsdErr);
	}
	return fail;
}

bool Mtrx::isIdentityMatrix(Double thresh) const
{
	bool fail = false;
	Double largestDiffOffdiagonal = ZERO;
	Double largestDiffDiagonal = ZERO;
	for(UInt j=0; j<getCol(); j++) {
		for(UInt i=0; i<getRow(); i++) {
			if (i != j && fabs(val(i,j)) > thresh) {
				fail = true;
				Double diff = fabs(val(i,j));
				if (diff>largestDiffOffdiagonal) largestDiffOffdiagonal = diff;
			}else if (i == j && fabs(val(i,j)-ONE) > thresh) {
				fail = true;
				Double diff = fabs(val(i,j)-ONE);
				if (diff>largestDiffDiagonal) largestDiffDiagonal = diff;
			}
		}
	}
	if (fail) {
		printf("the largest absolute difference for off-diagonal elements is %-16.10f\n",
				largestDiffOffdiagonal);
		printf("the largest absolute difference for diagonal elements is %-16.10f\n",largestDiffDiagonal);
		return false;
	}
	return true;
}
