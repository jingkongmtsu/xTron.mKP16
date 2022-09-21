/**
 * \file    molegrids.h
 * \brief   describing the grid system for molecule
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef MOLEGRIDS_H
#define MOLEGRIDS_H
#include <vector>
#include "libgen.h"
#include "atomgrids.h"
#include "scalablevec.h"

namespace xcintsinfor {
	class XCIntJobInfor;
}

namespace shellsize {
	class MolShellSize;
}

namespace shell {
	class MolShell;
}

namespace molegrids {

	using namespace atomgrids;
	using namespace xcintsinfor;
	using namespace shellsize;
	using namespace shell;
	using namespace std;

	///
	/// here we define the vector of atom grids
	///
	typedef std::vector<AtomGrids,tbb::scalable_allocator<AtomGrids> > MolGridVec;

	/**
	 * Aggregation of the atom grids will lead to the molecule grids
	 */
	class MoleGrids {

		private:

			MolGridVec molGrids;     ///< grid information 

			///
			/// the batchInfor stores the batch information for the given molecule
			/// the data fields simply divided into three fields: 
			///
			/// - the first one is the atom index of this batch;
			/// - the second one is the atomic number corresponding to this batch;
			/// - the third one is the local batch index for the AtomGrid
			///
			/// by using the above three information, we can build the batchGrid
			///
			UIntVec batchInfor;      

		public:

			/**
			 * constructor for building molecule grids in terms of atom types
			 */
			MoleGrids(const MolShell& ms, const MolShellSize& size, const XCIntJobInfor& infor);

			///
			/// destructor
			/// 
			~MoleGrids() { };

			/**
			 * number of terms in the molgrids array
			 */
			UInt getNAtomGrids() const {
				return molGrids.size();
			};

			/**
			 * get the maximum number of batch grid number among all of possible batch
			 */
			UInt getNMaxBatchGrids() const;

			/**
			 * get the total number of batches
			 */
			UInt getNTotalBatches() const { return batchInfor.size()/3; };

			/**
			 * for the given batch we return the atomic number as well as the 
			 * the local batch index
			 */
			void getBatchInfor(const UInt& globalBatchIndex, UInt& iAtomIndex, 
					UInt& Z, UInt& localBatchIndex) const {
				iAtomIndex      = batchInfor[3*globalBatchIndex  ];
				Z               = batchInfor[3*globalBatchIndex+1];      
				localBatchIndex = batchInfor[3*globalBatchIndex+2];      
			};

			/**
			 * get the corresponding atomGrids data from the given atomic number
			 */
			const AtomGrids& getAtomGrid(const UInt& Z) const;

			/**
			 * debug printing
			 */
			void print(UInt level = 1) const;

	};

}

#endif
