/**
 * \file   atomfockmtrx.h
 * \brief  this class stores the fock matrix pieces for J/K results in atom shell dimension
 * \author Fenglai Liu
 */
#ifndef ATOMFOCKMTRX_H
#define ATOMFOCKMTRX_H
#include "libgen.h"
#include "blockmatrix.h"
using namespace blockmatrix;

namespace shellpair {
	class AtomShellPair;
	class ShellPair;
}

namespace symmjkdigest {
	class SymmJKDigest;
}

namespace gintsinfor {
	class GIntJobInfor;
}

namespace blockmatrixlist {
	class BlockMtrxList;
}

namespace atomfockmtrx {

	using namespace shellpair;
	using namespace gintsinfor; 
	using namespace symmjkdigest;
	using namespace blockmatrixlist;

	///
	/// pre-defined the position macro for linking the atom shell result
	/// and the results for shell 
	///
	/// because for real shell quartet calculation shell pair may be switched,
	/// therefore we need to know how to link the result in atom shell and 
	/// the result in shell. For example, (12|34) may be finally calculated
	/// as (34|21), so the result k13 in terms of shell dimension (in symmJKDigest) 
	/// actually corresponding to atom shell result K23 in transpose form here 
	/// in this class. That's why we need to set up the linking information below.
	///
	const UInt LINKING_NULL = -1;
	const UInt LINKING_K13  =  1;  ///< the result in shell dimension corresponding to aK13/bK13 below
	const UInt LINKING_K14  =  2;  ///< the result in shell dimension corresponding to aK14/bK14 below
	const UInt LINKING_K23  =  3;  ///< the result in shell dimension corresponding to aK23/bK23 below
	const UInt LINKING_K24  =  4;  ///< the result in shell dimension corresponding to aK24/bK24 below
	const UInt LINKING_J12  =  5;  ///< the result in shell dimension corresponding to J12 below
	const UInt LINKING_J34  =  6;  ///< the result in shell dimension corresponding to J34 below


	///
	/// this is the class responsible to collect the Fock matrix results for JK digestion.
	///
	/// If the target system becomes very big, each working thread can not afford to 
	/// have to have its own full matrix result. In this case we can not use TBB_GInts4DMatrix,
	/// where each thread has its own full Fock matrix; instead we need to use the 
	/// class of TBB_GInts4DMtrxList, where each thread use block matrix list to store
	/// the fock matrix result.
	///
	/// this class is designed to work with TBB_GInts4DMtrxList. It stores the fock matrix
	/// result in terms of atom shell dimension after J/K digestion, and later it will be 
	/// pushed back to the block matrix list in TBB_GInts4DMtrxList.
	///
	class AtomFockMtrx {

		private:

			//
			// define the linking information for exchange results 
			// k13, k24 etc. are results stored in symmjkdigest class
			//
			bool k13InTrans;        ///< whether the content of k13 in transpose of the result matrix here?
			bool k14InTrans;        ///< whether the content of k14 in transpose of the result matrix here?
			bool k23InTrans;        ///< whether the content of k23 in transpose of the result matrix here?
			bool k24InTrans;        ///< whether the content of k24 in transpose of the result matrix here?
			UInt k13Pos;            ///< for result of k13 in shell dimension, who it corresponding to?
			UInt k14Pos;            ///< for result of k14 in shell dimension, who it corresponding to?
			UInt k23Pos;            ///< for result of k23 in shell dimension, who it corresponding to?
			UInt k24Pos;            ///< for result of k24 in shell dimension, who it corresponding to?

			//
			// define the linking information for coulomb results 
			// j12, j34 are results stored in symmjkdigest class
			//
			bool j12InTrans;        ///< whether the content of j12 in transpose of the result matrix here?
			bool j34InTrans;        ///< whether the content of j34 in transpose of the result matrix here?
			UInt j12Pos;            ///< for result of j12 in shell dimension, who it corresponding to?
			UInt j34Pos;            ///< for result of j34 in shell dimension, who it corresponding to?

			//
			// exchange result matrices 
			//
			BlockMtrx aK13;         ///< P24*(12|34) = K13 for alpha state
			BlockMtrx aK14;         ///< P23*(12|34) = K14 for alpha state
			BlockMtrx aK23;         ///< P14*(12|34) = K23 for alpha state
			BlockMtrx aK24;         ///< P13*(12|34) = K24 for alpha state
			BlockMtrx bK13;         ///< P24*(12|34) = K13 for beta  state
			BlockMtrx bK14;         ///< P23*(12|34) = K14 for beta  state
			BlockMtrx bK23;         ///< P14*(12|34) = K23 for beta  state
			BlockMtrx bK24;         ///< P13*(12|34) = K24 for beta  state

			//
			// coulomb result matrices
			//
			BlockMtrx J12;          ///< P34*(12|34) = J12
			BlockMtrx J34;          ///< P12*(12|34) = J34

		public:

			///
			/// constructor, doing allocation for memory for the matrix inside
			///
			AtomFockMtrx(UInt maxNBasis):k13InTrans(false),k14InTrans(false),
			k23InTrans(false),k24InTrans(false),
			k13Pos(LINKING_NULL),k14Pos(LINKING_NULL),k23Pos(LINKING_NULL),k24Pos(LINKING_NULL),
			j12InTrans(false),j34InTrans(false),j12Pos(LINKING_NULL),j34Pos(LINKING_NULL),
			aK13(maxNBasis,maxNBasis),aK14(maxNBasis,maxNBasis), 
			aK23(maxNBasis,maxNBasis),aK24(maxNBasis,maxNBasis),
			bK13(maxNBasis,maxNBasis),bK14(maxNBasis,maxNBasis),
			bK23(maxNBasis,maxNBasis),bK24(maxNBasis,maxNBasis),        
			J12 (maxNBasis,maxNBasis),J34 (maxNBasis,maxNBasis) { } 

			///
			/// destructor
			///
			~AtomFockMtrx() { };

			///
			/// initialize the information based on the atom input shell pair data
			/// and initialize the block form of density matrix data
			///
			void initInfor(const GIntJobInfor& infor, 
					const AtomShellPair& braAtomSP,const AtomShellPair& ketAtomSP);

			///
			/// reset the link information according to the shell pair information
			///
			void resetLinkInfor(bool switchSP, const ShellPair& braSP, const ShellPair& ketSP);

			///
			/// update the result from the local results
			///
			void updateResults(const GIntJobInfor& infor, const SymmJKDigest& localResults);

			///
			/// merge the results into the block matrix list
			/// iSpin identify the list spin state
			///
			void merge(const GIntJobInfor& infor, BlockMtrxList& list, UInt iSpin) const;
	};

};

#endif
