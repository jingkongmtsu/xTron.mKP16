/**
 * \file    dftmatrix.h
 * \brief   produce matrix result based on signifacant basis set order
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef DFTMATRIX_H
#define DFTMATRIX_H
#include "libgen.h"
#include "spinmatrix.h"
using namespace spinmatrix;

namespace sigatombasis {
	class SigAtomBasis;
}

namespace shell {
	class AtomShell;
	class MolShell;
}

namespace dftmatrix {

	//
	// define the work status of DFT matrix
	//
	const UInt ROW_IN_SIG_ORDER     = 1;
	const UInt COL_IN_SIG_ORDER     = 2;
	const UInt ROW_COL_IN_SIG_ORDER = 3;

	using namespace shell;
	using namespace sigatombasis;

	class DFTMatrix : public SpinMatrix {

		protected:

			UInt status;         ///< characterizes the status of DFT matrix (see above)

		public:

			/**
			 * constructor a partial DFT matrix - whose row is full size 
			 * of shell, column is significant order of basis set which 
			 * is characterized by colSigList
			 */
			DFTMatrix(const MolShell& rowShell, const SigAtomBasis& colSigList, const UInt& nSpin);

			/**
			 * constructor a partial DFT matrix - whose row is in significant 
			 * order of basis set which is characterized by rowSigList, and 
			 * column is full size of shell
			 */
			DFTMatrix( const SigAtomBasis& rowSigList, const MolShell& colShell, const UInt& nSpin);

			/**
			 * constructor a full DFT matrix - both of row and column are in
			 * significant order
			 */
			DFTMatrix( const SigAtomBasis& rowSigList, const SigAtomBasis& colSigList, const UInt& nSpin);

			/**
			 * constructor for building a "ATOM" based DFT matrix
			 * the row and column is nSigBasis for the given row Atom 
			 * and col Atom. Both of rowAtom and colAtom are local
			 * significant atom index
			 */
			DFTMatrix( const SigAtomBasis& rowSigList, const SigAtomBasis& colSigList,
					const UInt& rowAtom, const UInt& colAtom, const UInt& nSpin);

			/**
			 * destructor
			 */
			~DFTMatrix() { };

			/**
			 * this function is used to convert some "Atom" based AO matrix 
			 * into the significant order
			 * for example, convert the free atom density matrix into significant 
			 * order based on the given sigList information
			 * \param rowAtomShell, colAtomShell: for AO matrix, it's corresponding row 
			 *                                    and column atom shell data
			 * \param rowBasis,colBasis         : sigList information for both row and column
			 * \param iRowSigAtom,iColSigAtom   : local significant atom index for row and column
			 * \param M                         : input atom based AO matrix
			 *
			 */
			void intoSigOrder(const UInt& iSpin, const Mtrx& M, const AtomShell& rowAtomShell, 
					const AtomShell& colAtomShell, const SigAtomBasis& rowBasis, 
					const SigAtomBasis& colBasis, const UInt& iRowSigAtom,
					const UInt& iColSigAtom);

			/**
			 * This function is used to convert the matrix data into significant order
			 * this function utilizes two sig basis, therefore it could be used to convert
			 * both of the symmetrical/assymetrial matrix
			 * \param  iSpin : which spin state it refers?
			 * \param  source: matrix converting into sig order
			 * \param  b1    : sig basis set infor corresponding to the row conversion
			 * \param  b2    : sig basis set infor corresponding to the column conversion
			 */
			void intoSigOrder(const UInt& iSpin, const Mtrx& source, 
					const SigAtomBasis& b1, const SigAtomBasis& b2); 

			/**
			 * This function is used to convert the matrix data into significant order
			 * this function utilizes two sig basis, therefore it could be used to convert
			 * both of the symmetrical/assymetrial matrix
			 * \param  goal  : matrix converting into normal shell order
			 * \param  b1    : sig basis set infor corresponding to the row conversion
			 * \param  b2    : sig basis set infor corresponding to the column conversion
			 */
			void intoNormOrder(const SigAtomBasis& b1, const SigAtomBasis& b2, 
					const UInt& iSpin, Mtrx& goal) const;

			/**
			 * This function is used to convert the matrix data into significant order
			 * for row or column of the matrix
			 * \param  M     : matrix converting into sig order
			 * \param  b     : sig basis set information corresponding to the conversion
			 */
			void intoSigOrder(const UInt& iSpin, const Mtrx& M, const SigAtomBasis& b);

			/**
			 * This function is used to convert the matrix data into significant order
			 * for row or column of the matrix
			 * \param  M     : matrix converting into normal shell order
			 * \param  b     : sig basis set information corresponding to the conversion
			 */
			void intoNormOrder(const SigAtomBasis& b, const UInt& iSpin, Mtrx& M) const;

			/**
			 * debug purpose
			 */
			void print() const;
	};

	class DFTVect {

		protected:

			UInt nSpin;                ///< number of spin states
			UInt nData;                ///< data length
			DoubleVec vectList;        ///< data in terms of spin states, in vector form

		public:

			/**
			 * constructor a DFT vector
			 */
			DFTVect(const SigAtomBasis& sigList, const UInt& nSpin);

			///
			/// destructor
			///
			~DFTVect() { };

			/**
			 * This function is used to convert the vector into significant order
			 * \param  source: vector before converting
			 * \param  b     : sig basis set information corresponding to conversion
			 */
			void intoSigOrder(const DoubleVec& v, const SigAtomBasis& b);

			/**
			 * This function is used to convert the vector into normal order
			 * \param  source: vector before converting
			 * \param  b     : sig basis set information corresponding to conversion
			 */
			void intoNormOrder(DoubleVec& v, const SigAtomBasis& b) const;

			/**
			 * get the number of spin states
			 */
			UInt getNSpin() const { return nSpin; };

			/**
			 * return the vector data in const way
			 */
			const Double* getConstVect(const UInt& iSpin) const {
				UInt state = iSpin;
				if (nSpin == 1 && iSpin == 1) state = 0;
				return &vectList[state*nData];
			};

			/**
			 * debug purpose
			 */
			void print() const;
	};

}

#endif
