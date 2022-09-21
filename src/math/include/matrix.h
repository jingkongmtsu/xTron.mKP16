/**
 * \file   matrix.h
 * \author Fenglai Liu
 *
 * A wrapper function for matrix based on STL vector and mathematical
 * libraries (BLAS and LAPACK)
 *
 */
#ifndef MATRIX_H
#define MATRIX_H
#include<vector>
#include<string>
#include<boost/lexical_cast.hpp>
#include "libgen.h"
#include "excep.h"
#include "tbb/tbb.h"
#include "scalablevec.h"

namespace matrix {

	using namespace excep;
	using namespace tbb;
	using namespace std;

	//
	// declare the matrix class
	//
	class Mtrx;

	///
	/// utility functions for matrix
	///
	
	///
	/// if the input index is in pack form (that is to say,
	/// index goes from 0 to n*(1+n)/2); then we can use this
	/// function to decode the normal row/col index from the 
	/// pack form global index
	/// \param n    : row/column for the pack form data
	/// \param index: input global index in pack form
	/// \return rowIndex, colIndex: index corresponding to normal form of data
	///
	void getRowColIndexFromPackIndex(const UInt& n, const UInt& index, 
			UInt& rowIndex, UInt& colIndex);

	///
	/// here we define the vector of matrices
	///
	typedef std::vector<Mtrx,tbb::scalable_allocator<Mtrx> > MtrxVec;

	///
	/// class for normal matrix
	///
	class Mtrx {

		protected:

			UInt row;                ///< row dimension
			UInt col;                ///< column dimension
			UInt ld;                 ///< leading dimension
			DoubleVec array;         ///< used to hold the memory for data

		public:

			/****************************************************
			 *               CONSTRUCTORS etc.                  *
			 *also includes dimension change/memory manipulation*
			 ****************************************************/

			/**
			 * construct matrix with given row and column
			 */
			Mtrx(const UInt& r, const UInt& c):row(r),col(c),ld(r) {
				array.assign(row*col,ZERO);
			};

			/**
			 * initilizing matrix with given row and column
			 */
			void init(const UInt& r, const UInt& c) {
				row = r;
				col = c;
				ld  = r;
				array.assign(row*col,ZERO);
			};		

			/**
			 * default empty matrix class
			 */
			Mtrx():row(0),col(0),ld(0){  };

			/**
			 * copy constructor and assign operator
			 */
			Mtrx(const Mtrx& A);

			/**
			 * assign operator
			 */
			Mtrx& operator=(const Mtrx& A);

			///
			/// destructor
			///
			~Mtrx() {  }

			/**
			 * reset function is usually used when the matrix wants
			 * to reset it's dimension. This could lead to shrink
			 * or enlarge the current matrix. 
			 *
			 * \param nRow,nCol: the new dimension for the matrix
			 * \param reform: whether we will reforming the data
			 */
			void reset(const UInt& nRow, const UInt& nCol, bool reform);

			/**
			 * release the memory hold by the matrix
			 * this function should be called with CAUTION!!!
			 */
			void memRelease();

			/**
			 * enlarge the matrix to the given dimension
			 * \param r: the final matrix's row would be row+r
			 * \param c: the final matrix's col would be col+c
			 */
			void enlarge(const UInt& r, const UInt& c);

			//
			// clear the content of matrix
			// however, we may keep the dimension etc. information
			// since such data may need to be re-used
			//
			void clear() { array.assign(array.size(),ZERO); };

			/****************************************************
			 * inline operations, basically member data fetching*
			 ****************************************************/
			/**
			 * return the vector form of matrix in const way
			 */
			const DoubleVec& getVec() const { return array; };

			/**
			 * return the vector form of matrix in non-const way
			 */
			DoubleVec& getVec() { return array; };

			/**
			 * i is the index for row, j is the index for col
			 * such expression is used for a non-const matrix object
			 */
			Double& operator()(const UInt& i, const UInt& j) {
#ifdef DEBUG
				if((i>=row) || (j>=col)) {
					string info1 = "input row:" + boost::lexical_cast<string>(i); 
					string info2 = " input col:" + boost::lexical_cast<string>(j); 
					string info3 = " matrix row: " + boost::lexical_cast<string>(row);
					string info4 = " matrix col: " + boost::lexical_cast<string>(col);
					string infor = info1 + info2 + info3 + info4; 
					Excep excep("Mtrx","operator()",EXCEPTION_MATRIX_DIMENSION_OVERFLOW,infor);
					handleExcep(excep);
				}
#endif
				return array[i+j*ld];
			};

			/**
			 * i is the index for row, j is the index for col
			 * such expression is used for a const matrix object
			 */
			Double operator()(const UInt& i, const UInt& j) const {
#ifdef DEBUG
				if((i>=row) || (j>=col)) {
					string info1 = "input row:" + boost::lexical_cast<string>(i); 
					string info2 = " input col:" + boost::lexical_cast<string>(j); 
					string info3 = " matrix row: " + boost::lexical_cast<string>(row);
					string info4 = " matrix col: " + boost::lexical_cast<string>(col);
					string infor = info1 + info2 + info3 + info4; 
					Excep excep("Mtrx","operator(), const",EXCEPTION_MATRIX_DIMENSION_OVERFLOW,infor);
					handleExcep(excep);
				}
#endif
				return array[i+j*ld];
			};

			/**
			 * This function is actually same with the above operator()
			 * however, this one could be used when we need to refer
			 * this(i,j) situation.
			 */
			Double val(const UInt& i, const UInt& j) const {
#ifdef DEBUG
				if((i>=row) || (j>=col)) {
					string info1 = "input row:" + boost::lexical_cast<string>(i); 
					string info2 = " input col:" + boost::lexical_cast<string>(j); 
					string info3 = " matrix row: " + boost::lexical_cast<string>(row);
					string info4 = " matrix col: " + boost::lexical_cast<string>(col);
					string infor = info1 + info2 + info3 + info4; 
					Excep excep("Mtrx","val",EXCEPTION_MATRIX_DIMENSION_OVERFLOW,infor);
					handleExcep(excep);
				}
#endif
				return array[i+j*ld];
			};

			/**
			 * get the const pointer for the element in the array
			 * r is the position in row, c is the position in column
			 */
			const Double* getPtr(UInt r=0, UInt c=0) const {
#ifdef DEBUG
				if((r>=row) || (c>=col)) {
					string info1 = "input row:" + boost::lexical_cast<string>(r); 
					string info2 = " input col:" + boost::lexical_cast<string>(c); 
					string info3 = " matrix row: " + boost::lexical_cast<string>(row);
					string info4 = " matrix col: " + boost::lexical_cast<string>(col);
					string infor = info1 + info2 + info3 + info4; 
					Excep excep("Mtrx","getPtr, const",EXCEPTION_MATRIX_DIMENSION_OVERFLOW,infor);
					handleExcep(excep);
				}
#endif
				return &array[r+c*ld];
			};

			/**
			 * get the pointer for the element in the array
			 * r is the position in row, c is the position in column
			 */
			Double* getPtr(UInt r=0, UInt c=0) {
#ifdef DEBUG
				if((r>=row) || (c>=col)) {
					string info1 = "input row:" + boost::lexical_cast<string>(r); 
					string info2 = " input col:" + boost::lexical_cast<string>(c); 
					string info3 = " matrix row: " + boost::lexical_cast<string>(row);
					string info4 = " matrix col: " + boost::lexical_cast<string>(col);
					string infor = info1 + info2 + info3 + info4; 
					Excep excep("Mtrx","getPtr",EXCEPTION_MATRIX_DIMENSION_OVERFLOW,infor);
					handleExcep(excep);
				}
#endif
				return &array[r+c*ld];
			};

			///
			/// return the row of matrix
			///
			UInt getRow() const { return row; };

			///
			/// return the column of matrix
			///
			UInt getCol() const { return col; };

			///
			/// return the leading dimension of matrix
			///
			UInt getLd() const { return ld; };

			///
			/// whether the matrix is square?
			///
			bool isSquare() const { return (row == col); };

			///
			/// return the vector length
			///
			UInt vecLen() const { return array.size(); };

			/****************************************************
			 *             MATRIX COPY OPERATION                *
			 ****************************************************/

			/**
			 * this function is used to copy the matrix all to all
			 */
			void copyMatrix(const Mtrx& A);  

			/**
			 * copy data from the part of matrix to here
			 * \param rpos1:  row position for reading from the matrix(starting from 0)
			 * \param cpos1:  column position for reading from the matrix(starting from 0)
			 * \param rpos2:  row position for writing into matrix(starting from 0)
			 * \param cpos2:  column position for writing into matrix(starting from 0)
			 * \param  r   :  number of rows in matrix
			 * \param  c   :  number of columns in matrix
			 * \param  A   :  matrix for data source
			 */
			void copyFromMatrix(const UInt& rpos1, const UInt& cpos1, 
					const UInt& rpos2, const UInt& cpos2, 
					const UInt& r, const UInt& c, const Mtrx& A);

			/**
			 * copy data from the part of matrix to another matrix of A
			 * \param rpos1:  row position for reading from the matrix(starting from 0)
			 * \param cpos1:  column position for reading from the matrix(starting from 0)
			 * \param rpos2:  row position for writing into matrix(starting from 0)
			 * \param cpos2:  column position for writing into matrix(starting from 0)
			 * \param  r   :  number of rows in matrix
			 * \param  c   :  number of columns in matrix
			 * \param  A   :  matrix where we copy data to
			 */
			void copyToMatrix(const UInt& rpos1, const UInt& cpos1, 
					const UInt& rpos2, const UInt& cpos2, 
					const UInt& r, const UInt& c, Mtrx& A) const;

			/**
			 * update data from the part of matrix to here
			 * compared with the above copy function, here what we do is
			 * M(i,j) = M(i,j) + a*A(i,j)
			 * M is this matrix
			 * \param rpos1:  row position for reading from the matrix(starting from 0)
			 * \param cpos1:  column position for reading from the matrix(starting from 0)
			 * \param rpos2:  row position for writing into matrix(starting from 0)
			 * \param cpos2:  column position for writing into matrix(starting from 0)
			 * \param  r   :  number of rows in matrix
			 * \param  c   :  number of columns in matrix
			 * \param  A   :  matrix for data source
			 */
			void updateFromMatrix(const UInt& rpos1, const UInt& cpos1, 
					const UInt& rpos2, const UInt& cpos2, 
					const UInt& r, const UInt& c, const Double& a, const Mtrx& A);

			/**
			 * copy data from the part of matrix to another matrix of A
			 * compared with the above copy function, here what we do is
			 * A(i,j) = A(i,j) + a*M(i,j)
			 * M is this matrix
			 * \param rpos1:  row position for reading from the matrix(starting from 0)
			 * \param cpos1:  column position for reading from the matrix(starting from 0)
			 * \param rpos2:  row position for writing into matrix(starting from 0)
			 * \param cpos2:  column position for writing into matrix(starting from 0)
			 * \param  r   :  number of rows in matrix
			 * \param  c   :  number of columns in matrix
			 * \param  A   :  matrix where we copy data to
			 */
			void updateToMatrix(const UInt& rpos1, const UInt& cpos1, 
					const UInt& rpos2, const UInt& cpos2, 
					const UInt& r, const UInt& c, const Double& a, Mtrx& A) const;

			/**
			 * copy data from a vector to here for a sub-matrix
			 * \param rpos :  row position this matrix to recieve data(starting from 0)
			 * \param cpos :  column position for this matrix to recieve data(starting from 0)
			 * \param  r   :  number of rows in matrix
			 * \param  c   :  number of columns in matrix
			 * \param  v   :  the source vector data
			 */
			void copyFromVector(const UInt& rpos, const UInt& cpos, 
					const UInt& r, const UInt& c, const Double* v);

			/**
			 * copy data from the part of matrix to a vector
			 * \param rpos :  row position for reading from the matrix(starting from 0)
			 * \param cpos :  column position for reading from the matrix(starting from 0)
			 * \param  r   :  number of rows in matrix
			 * \param  c   :  number of columns in matrix
			 * \param  v   :  result vector where we copy data to
			 */
			void copyToVector(const UInt& rpos, const UInt& cpos, 
					const UInt& r, const UInt& c, Double* v) const;

			/**
			 * update data from a vector to a sub-matrix here
			 * compared with the above copy function, here what we do is
			 * M(i,j) = M(i,j) + a*v(k)
			 * M is this matrix
			 * \param rpos :  row position this matrix to recieve data(starting from 0)
			 * \param cpos :  column position for this matrix to recieve data(starting from 0)
			 * \param  r   :  number of rows in matrix
			 * \param  c   :  number of columns in matrix
			 * \param  v   :  vector for data source
			 */
			void updateFromVector(const UInt& rpos, const UInt& cpos, 
					const UInt& r, const UInt& c, const Double& a, const Double* v);

			/**
			 * copy data from the part of matrix to a vector
			 * compared with the above copy function, here what we do is
			 * v(k) = v(k) + a*M(i,j)
			 * M is this matrix
			 * \param rpos :  row position for reading from the matrix(starting from 0)
			 * \param cpos :  column position for reading from the matrix(starting from 0)
			 * \param  r   :  number of rows in matrix
			 * \param  c   :  number of columns in matrix
			 * \param  v   :  vector where we copy data to
			 */
			void updateToVector(const UInt& rpos, const UInt& cpos, 
					const UInt& r, const UInt& c, const Double& a, Double* v) const;

			/**
			 * if the given matrix is symmetrical, form the vector
			 * that contains its lower-triangular part plus the 
			 * diagonal part
			 */
			void formSymmVec(DoubleVec& vec) const;

			/**
			 * loading the matrix content from a given symmetrical
			 * vector
			 */
			void loadSymmVec(const DoubleVec& vec);

			/****************************************************
			 *        MATRIX OPERATION WITH ITSELF              *
			 ****************************************************/

			/**
			 * copy lower triangular part to upper triangular part
			 * this function should be for debugging purpose
			 */
			void copyLowToUpper();

			/**
			 * perform operation of A = A + A^{T}
			 */
			void addTranspose();

			/**
			 * do transpose operation via additional memory
			 */
			void simpleTranspose();

			/**
			 * do transpose operation in place
			 */
			void transposeInPlace();

			/**
			 * interface function for transposing matrix
			 */
			void transpose(bool withInPlaceTranspose) {
				if (withInPlaceTranspose) {
					transposeInPlace();
				}else{
					simpleTranspose();
				}
			};

			/****************************************************
			 *             BASIC MATRIX OPERATIONS - BLAS1      *
			 *         here the operations are parallelized     *
			 ****************************************************/
			/**
			 * A = a*A
			 */
			void scale(Double a);

			/**
			 * A = a*M + A
			 */
			void add(const Mtrx& M, Double a=ONE);

			/**
			 * set the whole matrix to some certain value of a
			 */
			void set(Double a);

			/**
			 * this is to calculate the rmsd difference between
			 * this matrix and the given matrix A
			 *
			 * if onlyLowerPart is passed in, then we only use
			 * A and this matrix's lower part
			 */
			Double rmsd(const Mtrx& A, bool onlyLowerPart) const;

			/****************************************************
			 *      BLAS TYPE MATRIX OPERATIONS -BLAS 2/3       *
			 *       here the operations are parallelized       *
			 ****************************************************/
			/**
			 * matrix multiplication
			 * C = alpha*A*B + beta*C or C = alpha*B*A + beta*C
			 * A is symmetrical matrix, so only it's lower triangular part is referred
			 * if A*B, then state is L
			 * else it's B*A state is R
			 */
			void symMatrixMult(const Mtrx& A, const Mtrx& B, const char& state,
					Double alpha, Double beta, bool withPara = false);

			/**
			 * matrix multiplication
			 * we do C = alpha*C*A + beta*C or C = alpha*A*C + beta*C
			 * here inside we will judge whether the matrix A is symmetrical
			 * if A is symmetrical, then we will call msymmul
			 * else we will call mmul
			 * \param isSymm: wheather the A is symmetrical
			 * \param state: if C*A then A is on right(R), else A is on left(L)
			 */
			void sideMult(const Mtrx& A, const char& state, bool isSymm,
					Double alpha, Double beta, bool withPara = false); 

			/**
			 * matrix multiplication for the full matrix. 
			 * this is for non-symmetrical case
			 * Here we do not take care of the transpose situation, it will be handled
			 * in the mmul function
			 */
			void mult(const Mtrx& A, const Mtrx& B, bool transposeA, 
					bool transposeB, Double alpha, Double beta, bool withPara = false);

			/**
			 * multiplication for part of matrix in A and B
			 * \param  A, B: the full matrix of A and B, can not be in pack form
			 * \param  rowAPos, colAPos, rowBPos, colBPos: positions to locate the sub matrix
			 * \param  rowA, colA, rowB, colB: dimensions for sub matrix
			 * we note that the result matrix could be in pack form or not
			 */
			void mult(const Mtrx& A, const Mtrx& B, const UInt& rowAPos, 
					const UInt& colAPos, const UInt& rowBPos, const UInt& colBPos,
					const UInt& rowA, const UInt& colA, const UInt& rowB, const UInt& colB, 
					bool transposeA, bool transposeB, Double alpha, 
					Double beta, bool withPara = false);

			/**
			 * wrapper for the dgemv
			 * y = \alpha op(A) x + \beta y
			 * \param withLowerTriAngularPart whether we perform work with only lower part
			 */
			void mmulv(const Double* x, Double* y, const Double& alpha, 
					const Double& beta, bool withLowerTriAngularPart, 
					bool needTranspose, bool withPara = false) const;

			/**
			 * doing direct summation for the given matrix
			 * \param Ai the sub matrix which is going to be inserted in the matrix
			 * \param rowPos the row position for inserting
			 * \param colPos the col position for inserting
			 */
			void addPieceInDirectSum(const Mtrx& Ai, const UInt& rowPos, 
					const UInt& colPos);

			/**
			 * this functions performs the result that 
			 * e = sum_{ij}A_{ij}D_{ij}
			 *
			 * depending on that whether the given matrix is symmetrical or 
			 * not, we will use different scheme
			 */
			Double dotProduct(const Mtrx& D, bool isSymm) const;

			/****************************************************
			 *         LAPACK TYPE MATRIX OPERATIONS            *
			 ****************************************************/
			/**
			 * inverse of the squre matrix in general type
			 * the content is going to be destroyed
			 */
			void inverseGenMatrix(bool withPara = false);

			/**
			 * inverse of the symmetrical positively defined matrix
			 * the content is going to be destroyed
			 */
			void inverseSymmMatrix(bool withPara = false); 

			/**
			 * get the matrix power of p, given matrix should be 
			 * symmetrical
			 *
			 * we will use the same method as shown in formOrthoMatrix
			 * to generate the matrix power
			 */
			void power(Double p, bool withPara, Double thresh);

			/**
			 * calculate the orthogonal matrix for the original symmetrical matrix
			 * according to the original overlap matrix ov, and linear dependency threshold
			 */
			void formOrthoMatrix(const Mtrx& ov, Double linearDepThresh, bool withPara = false); 

			/**
			 * get eigenvectors and eigen values for a symmetrical matrix 
			 * the eigenvectors will be stored into this matrix, thus
			 * the original content is destroyed
			 */
			void doEigen(DoubleVec& eigenValue, bool withPara = false);

			/**
			 * compute the matrix's rank to detect linear dependency
			 */
			UInt rank(Double linearDepThresh, bool withPara = false) const;

			/****************************************************
			 *                   OTHERS                         *
			 ****************************************************/

			/**
			 * print the matrix with title and pre-defined columns
			 */
			void print(string title, UInt cols = 6) const;

			/**
			 * this function is used to compare the elements between
			 * this matrix and input matrix A for a given threshold
			 * thresh
			 *
			 * onlyCompareLowPart indicates whether we compare only
			 * the lower part of matrix
			 *
			 * if printMore is true, then we will print out all of 
			 * difference value when it's larger than the given 
			 * input threshold.
			 */
			bool diffComp(const Mtrx& A, const Double& thresh, 
					bool onlyCompareLowPart, bool printMore=false) const;

			/**
			 * this function examples that whether this matrix 
			 * is identity matrix
			 *
			 * for debug purpose
			 */
			bool isIdentityMatrix(Double thresh) const;
	};

}

#endif

