/**
 * \file   blockmatrix.h
 * \author Fenglai Liu
 */
#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H
#include "matrix.h"
using namespace matrix;

namespace blockmatrix {

	///
	/// class for block matrix (submatrix)
	///
	/// row/col Pos represent the position for this submatrix in terms of it's mother matrix.
	/// This is the relation between this matrix and its mother matrix. Such relation should 
	/// be defined initially (in constructor or init function), and never changed even though
	/// block matrix content is changed. 
	///
	///
	class BlockMtrx : public Mtrx {

		protected:

			UInt rowPos;         ///< row position of this submatrix to it's mother, from 0
			UInt colPos;         ///< col position of this submatrix to it's mother, from 0

		public:

			/****************************************************
			 *               CONSTRUCTORS etc.                  *
			 ****************************************************/

			/**
			 * construct matrix with given row and column dimension
			 * It's used to reserve mem space for the block matrix object
			 *
			 * in this constructor, rowPos and colPos will be left to 
			 * be initialized later
			 */
			BlockMtrx(const UInt& r, const UInt& c):Mtrx(r,c),rowPos(-1),
			colPos(-1) { };

			/**
			 * construct matrix with all of information given 
			 * \param r     row of block matrix
			 * \param c     col of block matrix
			 * \param rPos  row position for block matrix in it's mother matrix
			 * \param cPos  col position for block matrix in it's mother matrix
			 */
			BlockMtrx(const UInt& r, const UInt& c, const UInt& rPos, 
					const UInt& cPos):Mtrx(r,c),rowPos(rPos),colPos(cPos) { };

			/**
			 * reset the block matrix data
			 */
			void resetBlockData() {
				clear();
				rowPos = -1;
				colPos = -1;
			};

			/**
			 * initilize matrix with all of information given 
			 *
			 * \param r     row of block matrix
			 * \param c     col of block matrix
			 * \param rPos  row position for block matrix in it's mother matrix
			 * \param cPos  col position for block matrix in it's mother matrix
			 */
			void init(const UInt& r, const UInt& c, const UInt& rPos, 
					const UInt& cPos) {
				rowPos       = rPos;
				colPos       = cPos;
				bool needReform = true;
				reset(r,c,needReform);
			};	

			/**
			 * copy constructor and assign operator
			 */
			BlockMtrx(const BlockMtrx& A) {
				row    = A.getRow();      
				rowPos = A.getRowPos();      
				col    = A.getCol();
				colPos = A.getColPos();
				// here we can not copy the leading dimension information
				// since we are unable to access to the ld of A
				ld     = row; 
				array.assign(row*col,ZERO);
				copyMatrix(A);
			};

			/**
			 * assign operator
			 */
			BlockMtrx& operator=(const BlockMtrx& A) {
				if (this == &A) {
					return *this;
				}else{
					row    = A.getRow();      
					rowPos = A.getRowPos();      
					col    = A.getCol();
					colPos = A.getColPos();
					// here we can not copy the leading dimension information
					// since we are unable to access to the ld of A
					ld     = row; 
					array.assign(row*col,ZERO);
					copyMatrix(A);
					return *this;
				}
			};

			/**
			 * bool operator == , testing that whether the input
			 * block matrix is same with this one
			 *
			 * we note that if the input matrix is in transpose form
			 * of this matrix, we still consider that they are same
			 */
			bool operator==(const BlockMtrx& M) const {

				// test that whether two matrix are identical
				if ((rowPos == M.getRowPos() && colPos == M.getColPos()) && 
						(row == M.getRow() && col == M.getCol())) { 
					return true;
				}

				// test that whether two matrix are same but one is 
				// transpose of the other
				if ((rowPos == M.getColPos() && colPos == M.getRowPos()) && 
						(row == M.getCol() && col == M.getRow())) { 
					return true;
				}

				// now two matrix are not same
				return false;
			};

			///
			/// destructor
			///
			~BlockMtrx() {  };


			/****************************************************
			 * inline operations, basically member data fetching*
			 ****************************************************/

			///
			/// transpose the matrix content as required
			/// also the row position and column position are changed accordingly
			///
			void transposeBlockMatrix(bool inPlace = true) {
				transpose(inPlace);
				UInt pos = rowPos;
				rowPos   = colPos;
				colPos   = pos;
			};

			///
			/// return the row position for this submatrix in its mother matrix
			///
			UInt getRowPos() const { return rowPos; };

			///
			/// return the col position for this submatrix in its mother matrix
			///
			UInt getColPos() const { return colPos; };

			///
			/// whether the block matrix from a diagonal part?
			///
			bool isDiagonal() const { return (rowPos == colPos); };

			///
			/// whether the block matrix is from an upper triangular part?
			///
			bool fromUpper() const {
				if (isDiagonal()) return false;
				return (rowPos<colPos);
			};

			///
			/// whether the block matrix is from an lower triangular part?
			///
			bool fromLower() const {
				if (isDiagonal()) return false;
				return (rowPos>colPos);
			};

			///
			/// whether both of two blocks are exactly same
			///
			bool same(const BlockMtrx& M) const {
				if ((rowPos == M.getRowPos() && colPos == M.getColPos()) && 
						(row == M.getRow() && col == M.getCol())) { 
					return true;
				}
				return false;
			};

			///
			/// whether both of two blocks are same in transpose status
			///
			bool sameInTranspose(const BlockMtrx& M) const {
				if ((rowPos == M.getColPos() && colPos == M.getRowPos()) && 
						(row == M.getCol() && col == M.getRow())) { 
					return true;
				}
				return false;
			};

			/****************************************************
			 *                 MATRIX OPERATION                 *
			 ****************************************************/
			/**
			 * get data from mother matrix to here
			 * \param  A   :  mother matrix
			 * \param  inTranspose: whether this block matrix is in 
			 *         transpose status comparing against it's mother
			 */
			void getData(const Mtrx& A, bool inTranspose);

			/**
			 * update data from this matrix to its mother matrix
			 * \param  A   :  mother matrix
			 * \param  inTranspose: whether this block matrix is in 
			 *         transpose status comparing against it's mother
			 */
			void updateData(Mtrx& A, bool inTranspose) const;

			//
			// adding data from the given block matrix
			// they should be the same block matrix then
			// it's able to perform adding
			//
			void add(const BlockMtrx& A);

			//
			// print the matrix
			// since we still want to call the print function in
			// matrix, this function here called blockPrint
			//
			void blockPrint(const string& title) const;
	};

}

#endif

