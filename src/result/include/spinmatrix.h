/**
 * \file    spinmatrix.h
 * \brief   define the matrix results which is generated in terms of spin states
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef SPINMATRIX_H
#define SPINMATRIX_H
#include "libgen.h"
#include "matrix.h"
using namespace matrix;

namespace spinmatrix {

	/**
	 * \class   SpinMatrix
	 * \brief   describe matrices in terms of spin states
	 *
	 * all of stuff is here in this head file
	 */
	class SpinMatrix {

		protected:

			MtrxVec  mtrxList;   ///< the result matrices in term of spins

		public:

			///
			/// constructor to build the blank matrices
			/// 
			SpinMatrix(const UInt& nSpins, const UInt& row, const UInt& col):mtrxList(nSpins) {
				for(UInt iSpin=0; iSpin<nSpins; iSpin++) {
					mtrxList[iSpin].init(row,col);
				}
			};

			///
			/// we add in default constructor
			/// we take care of the spin state, 
			/// but the data length of each matrix is zero
			///
			/// the data will be initilized later
			///
			SpinMatrix(UInt nSpins):mtrxList(nSpins) { };

			///
			/// init function to initilize the row and column of matrix
			///
			void init(const UInt& row, const UInt& col) {
				for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
					mtrxList[iSpin].init(row,col);
				}
			};

			///
			/// default destructor
			///
			~SpinMatrix() { };

			///
			/// return the number of spin states
			///
			UInt getNSpin() const { return mtrxList.size(); };

			///
			/// return the mtrx in const state
			///
			const Mtrx& getMtrx(const UInt& iSpin) const {
				UInt n = iSpin;
				if (n == 1 && getNSpin() == 1) n = 0;
				return mtrxList[n];
			};

			///
			/// return the mtrx in non-const state
			///
			Mtrx& getMtrx(const UInt& iSpin) {
				UInt n = iSpin;
				if (n == 1 && getNSpin() == 1) n = 0;
				return mtrxList[n];
			};

			///
			/// initilize the data with the value given 
			///
			void initialize(const Double& a) {
				for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
					mtrxList[iSpin].set(a);
				}
			};

			///
			/// print in debug purpose
			///
			void print(const string& title) const {
				for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
					mtrxList[iSpin].print(title);
				}
			};

			///
			/// whether the data matrix here is square?
			///
			bool isSquare() const {
				if (getNSpin() == 2) {
					return (mtrxList[0].isSquare() && mtrxList[1].isSquare());
				}else{
					return mtrxList[0].isSquare();
				}
			};
	};

}

#endif
