/**
 * \file    atomgrids.h
 * \brief   describing the atom grid used in the XC part
 * \author  Fenglai Liu and Jing Kong
 * \note 
 * 
 * SG1 grid reference:
 *
 * Gill, Peter MW, Benny G. Johnson, and John A. Pople. 
 * "A standard grid for density functional calculations." 
 * Chemical Physics Letters 209.5 (1993): 506-512.
 *
 * Baker grid is contributed from John Baker, from Peter Pulay's group.
 * This part of code is not published.
 */
#ifndef ATOMGRIDS_H
#define ATOMGRIDS_H
#include "libgen.h"
#include "scalablevec.h"
#include "gridinfor.h"
#include "tbb/tbb.h"
#define lebedev  lebedev_

/**
 * C wrapper function for making the lebedev grids
 * \param NPT:    number of points
 * \return X,Y,Z: angular point coordinates
 * \return W:     weights for angular point 
 */
extern "C" {
	void lebedev(Double* X, Double* Y, Double* Z, Double* W, Int* NPT);
}

namespace xcintsinfor {
	class XCIntJobInfor;
}

namespace atomgrids {

	using namespace xcintsinfor;
	using namespace gridinfor;

	//
	// declare of class
	//
	class BatchGridInfor;
	class AtomGrids;

	///
	/// here we define the vector of batch grid infor
	///
	typedef std::vector<BatchGridInfor,tbb::scalable_allocator<BatchGridInfor> > BatchGridInforVec;

	//////////////////////////////////////////////////////////////
	//      class used to hold the batch grid information       //
	//////////////////////////////////////////////////////////////
	/**
	 * This class is just to hold the batch grid information for 
	 * everything batch in atom grid
	 * It's only used in atom grid class therefore everything is public
	 */
	class BatchGridInfor {

		public:

			// batch information
			UInt batchSize;             ///< the size for current batch

			// radial points information
			UInt initRadIndex;          ///< the starting radial point index
			UInt nRadPts;               ///< how many radial points included in this batch

			/**
			 * for each radial point, it's angular grid offset
			 * for example, if the given radial point is in region 3,
			 * then angPts[3*angGridOffset[iRad]] will gives the 
			 * proper pointer for angular point array
			 * angWts[angGridOffset[iRad]] gives the pointer for 
			 * angular weights array
			 */
			UIntVec angGridOffset; 

			/**
			 * for each radial point, it has a angular point region
			 * and for this region, here it records  
			 * it's begin and end index. For each radial point,
			 * it could cover all of its angular points, and it 
			 * could cover part of it.
			 */
			UIntVec angGridIndex;  

			/**
			 * whether the last radial grid point has its full
			 * angular points. If yes, then next batch will
			 * continue to the next radial point, else the next
			 * batch will still stick to the last radial point
			 */
			bool lastPointAngGridInParts; 

			///
			/// constructor
			/// \param atomGrids current atom grid object
			/// \param iBatch    this batch index
			/// \param initBS    the initial batch size
			///
			BatchGridInfor(const AtomGrids& atomGrids, const UInt& iBatch, const UInt& initBS);

			///
			/// destructor
			///
			~BatchGridInfor() { };

			///
			/// get the radial point starting index through previous batch
			///
			UInt getRadBegin(const AtomGrids& atomGrids, const UInt& iBatch) const; 

			///
			/// get the angular point starting index through previous batch
			///
			UInt getAngBegin(const AtomGrids& atomGrids, const UInt& iBatch) const; 
	};

	/**
	 * \class   AtomGrids
	 * \brief   generating grid information for the atom
	 */
	class AtomGrids {

		private:

			// prune grid information
			UInt pruneMethod;        ///< the prune method 
			UInt choice;             ///< grid choice - standrad, coarse or fine?
			UInt nRegions;           ///< number of prune regions
			DoubleVec scaleFac;      ///< scale factor array for pruning 
			UIntVec angInRegion;     ///< the angular grid number for each region

			// general information for atom and shell
			UInt Z;                  ///< atomic number
			Double atomSize;         ///< atom size from shell data
			Double atomRadii;        ///< the radii for atom

			// general grid information
			UInt nTotalPts;          ///< the total number of points
			UInt nRadPts;            ///< number of radial points
			UInt nAngPts;            ///< number of angular points
			DoubleVec angPts;        ///< angular points coordinates 
			DoubleVec angWts;        ///< weights for the angular points
			DoubleVec radPts;        ///< radial points coordinates
			DoubleVec radWts;        ///< weights for the radial points

			// batch grid information
			UInt initBatchSize;      ///< input batch size for batch division
			BatchGridInforVec infor; ///< data to record the batch information  

			/////////////////////////////
			//     function calls      //
			/////////////////////////////

			/**
			 * check the input for making data
			 */
			void preCheck() const;

			/**
			 * form the grid information
			 */
			void formGridInfor();

			/**
			 * create the quadrature points and weights
			 */
			void formQuadraturePtsWts(const XCIntJobInfor& infor);

			/**
			 * create the batch
			 */
			void formBatch(const XCIntJobInfor& infor);

		public:

			/**
			 * constructor of the atom grid
			 * Generally, for a given atom(determined by the atomic number), as well
			 * as it's shell size(characterized by atom size); we are able to form
			 * the grid data for this atom.
			 */
			AtomGrids(const UInt& atomic, const Double& atomSize0, const XCIntJobInfor& xcInfor); 

			///
			/// destructor
			///
			~AtomGrids() { };

			///
			/// print for debug purpose
			///
			void print(UInt level = 0) const;


			///////////////////////////////////////////////////////
			// !!! functions related to the grid information     //
			///////////////////////////////////////////////////////

			/**
			 * whether we use pruned grids?
			 */
			bool usePruneGrids() const { 
				if (pruneMethod == NON_PRUNE_GRID) return false; 
				return true;
			};

			/**
			 * selecting the angular region by given the radius value
			 */
			UInt getRegion(const Double& radius) const;

			/**
			 *  get the radius value for the given rad index
			 */
			Double getRadius(const UInt& iRad) const { return radPts[iRad]; };

			/**
			 * return the number of ang points for the given region
			 */
			UInt getNAngForRegion(const UInt& region) const { 
				return angInRegion[region];
			};

			/**
			 * return the total number of angular points in terms 
			 * of angular grid type
			 */
			UInt getTotalAngPtsInType() const;

			/**
			 * get the total number of radial points
			 */
			UInt getNRadPts() const { return nRadPts; };

			/**
			 * get the total number of points
			 */
			UInt getNGrid() const { return nTotalPts; };

			/**
			 * get the offset for the given region in the angPts and angWts
			 * array
			 */
			UInt getAngGridOffset(const UInt& region) const;

			/**
			 * creating the radial grids by using Euler–Maclaurin formula
			 * see paper of 
			 * A standard grid for density functional calculations
			 * PMW Gill, BG Johnson, JA Pople - Chemical Physics Letters, 1993
			 *
			 * \param  nRad  number of radial points
			 * \return Pts   coordinate of the radial points
			 * \return Wts   weights for the radial points
			 */
			void eulmac(const UInt& nRad, Double* Pts, Double* Wts) const; 

			///////////////////////////////////////////////////////
			// !!! functions related to the atom information     //
			///////////////////////////////////////////////////////
			/**
			 * get the atomic number
			 */
			UInt getAtomic() const { return Z; };

			/**
			 * get atom size
			 */
			Double getAtomSize() const { return atomSize; };

			/**
			 * get atom radii
			 */
			Double getAtomRadii() const { return atomRadii; };

			///////////////////////////////////////////////////////
			// !!! functions related to the batch information    //
			///////////////////////////////////////////////////////
			/**
			 * number of batches 
			 */
			UInt getNBatch() const {
				return infor.size();
			};

			/**
			 * get the input batch size
			 * however, the batch size is not constant when we really
			 * divide the data into batches
			 * therefore the final batch size is given in batchGridInfor
			 * class
			 */
			UInt getInputBatchSize() const { return initBatchSize; };

			/**
			 * get the batch size for the given batch number
			 * This is the real batch size
			 */
			UInt getBatchSize(const UInt& iBatch) const { 
				return infor[iBatch].batchSize; 
			};

			/**
			 * return the current batch grid information 
			 */
			const BatchGridInfor& getBatchGridInfor(const UInt& iBatch) const{
				return infor[iBatch];
			};

			/**
			 * get the angular points for the given batch 
			 * and the given radial point
			 * iRad should be the relative index for this batch
			 * from 0 to nRadPts in bginfor
			 */
			const Double* getAngPts(const UInt& iBatch, const UInt& iRad) const {
				const BatchGridInfor& bgInfor = infor[iBatch];
				UInt offset = bgInfor.angGridOffset[iRad];
				return &angPts[3*offset];
			};

			/**
			 * get the angular grid index
			 * iRad should be the relative index for this batch
			 * from 0 to nRadPts in bginfor
			 */
			void getAngGridIndex(const UInt& iBatch, const UInt& iRad,
					UInt& begin, UInt& end) const {
				const BatchGridInfor& bgInfor = infor[iBatch];
				begin = bgInfor.angGridIndex[2*iRad];
				end   = bgInfor.angGridIndex[2*iRad+1];
			};

			/**
			 * get the angular point weights, similar to the getAngPts
			 * iRad should be the relative index for this batch
			 * from 0 to nRadPts in bginfor
			 */
			const Double* getAngWts(const UInt& iBatch, const UInt& iRad) const {
				const BatchGridInfor& bgInfor = infor[iBatch];
				UInt offset = bgInfor.angGridOffset[iRad];
				return &angWts[offset];
			};

			/**
			 * get the radial points for the given batch
			 */
			const Double* getRadPts(const UInt& iBatch) const {
				const BatchGridInfor& bgInfor = infor[iBatch];
				UInt radStart = bgInfor.initRadIndex;
				return &radPts[radStart];
			};

			/**
			 * get the radial point weights for the given batch
			 */
			const Double* getRadWts(const UInt& iBatch) const {
				const BatchGridInfor& bgInfor = infor[iBatch];
				UInt radStart = bgInfor.initRadIndex;
				return &radWts[radStart];
			};

			/**
			 * get the number of radial points in this batch
			 */
			UInt getNRadInBatch(const UInt& iBatch) const {
				const BatchGridInfor& bgInfor = infor[iBatch];
				return bgInfor.nRadPts;
			};
	};

}

#endif
