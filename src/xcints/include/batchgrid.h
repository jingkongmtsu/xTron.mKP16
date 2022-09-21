/**
 * \file    batchgrid.h
 * \brief   describing the batch grid 
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef BATCHGRID_H
#define BATCHGRID_H
#include "libgen.h"
#include "scalablevec.h"

namespace atom {
	class Atom;
}

namespace shell {
	class MolShell;
}

namespace sigatombasis {
	class SigAtomBasis;
}

namespace atomgrids {
	class AtomGrids;
}

namespace xcintsinfor {
	class XCIntJobInfor;
}

namespace batchgrid {

	using namespace atom;
	using namespace shell;
	using namespace sigatombasis;
	using namespace xcintsinfor;
	using namespace atomgrids;
	using namespace std;

	/**
	 * calculate the batch grids
	 */
	class BatchGrid {

		private:

			// grid data
			// the grid data is actually residing on a hollow disk
			// which centering on the batchAtom(parent atom). This
			// disk has a inside radius and outside radius
			Double inRad;           ///< inside radius
			Double outRad;          ///< outside radius
			UInt nGrids;            ///< number of grid points in this batch
			UInt batchAtom;         ///< for this batch of data, which parent atom it is
			DoubleVec coord;        ///< coordinates for batch grid points (ngrids*3)
			DoubleVec wts;          ///< weights for each grid point
			void init1Point(const AtomGrids& atomG, const Atom& atom, const UInt& iBatch);
			void init(const AtomGrids& atomG, const Atom& atom, const UInt& iBatch);
		public:

			///
			/// constructor for build the batch grid
			/// \param atomG0  atom grids data
			/// \param atom    provide the center's coordinates
			/// \param iBatch  the batch idex 
			/// 
			BatchGrid(const AtomGrids& atomG0, const Atom& atom, const UInt& iBatch, 
				const bool& onePoint = false);

			///
			/// destructor
			///
			~BatchGrid() { };

			/**
			 * print in debug purpose
			 */
			void print() const;

			/**
			 * return the batch grids
			 */
			UInt getNGrids() const {
				return nGrids;
			};

			/**
			 * return the coordinate values
			 */
			const Double* getGridCoord() const {
				return &coord[0];
			};

			/**
			 * return the weight values
			 */
			const Double* getGridWts() const {
				return &wts[0];
			};

			/**
			 * return the atom index in current batch
			 */
			UInt AtomInCurrentBatch() const {
				return batchAtom;
			};

			/**
			 * return the inside radius and outside radius
			 */
			void getBatchGridRadius(Double& inR, Double& outR) const {
				inR  = inRad;
				outR = outRad;
			};

			/******************************************************
			 *              Becke weights calculation             *
			 ******************************************************/
			/**
			 * this function is used to create the zero order Becke weights
			 * by using the atom coordinates in mol shell and the neighbor
			 * atom information from sigList
			 */
			void BeckeWeights0(const XCIntJobInfor& xcinfor, const MolShell& ms, 
					const SigAtomBasis& sigList); 
	};

}

#endif
