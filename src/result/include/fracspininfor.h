/**
 * \file    fractional spin calculation file
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef FRACSPININFOR_H
#define FRACSPININFOR_H
#include "libgen.h"
#include "molecule.h"
using namespace molecule;

namespace globalinfor {
	class GlobalInfor;
}

namespace fracspininfor {

	using namespace globalinfor; 

	/**
	 * this class collecting information for fractional spin calculation.
	 */
	class  FracSpinInfor{

		private:

			//
			// scaling mo relating information
			//
			UInt nAlphaScaledMO;        ///< the number of scaled mo for beta
			UInt alphaMOBeginIndex;     ///< from which alpha MO we begin to scale 
			DoubleVec alphaMOScalVals;  ///< values in a vector to scale alpha mo
			UInt nBetaScaledMO;         ///< the number of scaled mo for beta
			UInt betaMOBeginIndex;      ///< from which beta MO we begin to scale 
			Double betaMOScalVal;       ///< value to scale beta mo
			DoubleVec betaMOScalVals;   ///< value to scale beta mo
			Molecule newMol;            ///< new molecule information

		public:

			/**
			 * constructor
			 * \param  mol  molecule data
			 */
			FracSpinInfor(const GlobalInfor& globInfor,const Molecule& mol);

			/**
			 * default destructor
			 */
			~FracSpinInfor() { };

			/**
			 * do we do any fractional spin scaling job?
			 */
			bool doFracSpinJob() const { return (nAlphaScaledMO > 0 || nBetaScaledMO > 0); };

			/**
			 * return the modified molecule
			 */
			const Molecule& getMol() const { return newMol; };

			/**
			 * return the svaled value
			 */
			Double getScaleVal(UInt iSpin, UInt iMO) const { 
				if (iSpin == 0) return alphaMOScalVals[iMO]; 
				return betaMOScalVals[iMO]; 
			};

			/**
			 * return number of scaled mo
			 */
			UInt getNScaledMO(UInt iSpin) const { 
				if (iSpin == 0) return nAlphaScaledMO;   
				return nBetaScaledMO; 
			};

			/**
			 * get the beginning index of MO for scaling
			 */
			UInt getMOBeginIndex(UInt iSpin) const { 
				if (iSpin == 0) return alphaMOBeginIndex;  
				return betaMOBeginIndex; 
			}; 
	};
}

#endif

