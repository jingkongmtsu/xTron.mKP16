/**
 * \file    xcenergyinfor.h
 * \brief   this is used to analyze the XC energy component
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef XCENERGYINFOR_H
#define XCENERGYINFOR_H
#include <string>
#include "libgen.h"
#include "scalablevec.h"

namespace batchgrid {
	class BatchGrid;
}

namespace xcintsinfor {
	class XCIntsinfor;
}

namespace xcfunc {
	class XCFunc;
}

namespace xcenergyinfor {

	using namespace xcfunc;
	using namespace batchgrid;
	using namespace xcintsinfor;

	///
	/// here we define the component functional names that 
	/// used for energy profiling purpose
	///
	/// the reason we use the integer rather than string,
	/// is that we avoid to use data like 
	/// vector<string> xcENameList. Because this class
	/// will be used in multi-threading env, vector<string>
	/// is not though the default TBB allocator. Therefore
	/// we choose a more cumbersome way to avoid it.
	///
	const UInt XC_E_PROFILE_B13COOR_OPP        = 1;
	const UInt XC_E_PROFILE_B13COOR_PAR        = 2;
	const UInt XC_E_PROFILE_KP14C_U_STATIC_OPP = 3;
	const UInt XC_E_PROFILE_KP14C_U_STATIC_PAR = 4;

	/**
	 * \class   XCEnergyInfor
	 *
	 * this class is used to store the XC energy components, used for parameter optimization purpose etc.
	 */
	class XCEnergyInfor {

		private:

			UInt nGrids;           ///< the current batch number of grid points
			UInt nEComponents;     ///< number of XC energy components
			DoubleVec xcEList;     ///< energy components for the xc energy, accumulated from each batch
			UIntVec xcENameList;   ///< each energy component name, used macros defined above
			DoubleVec batchEList;  ///< stores the energy per grid on the current batch

			/**
			 * for the given functional name, return the position in xcENameList
			 */
			UInt getTheNameCode(const string& xcFuncName, UInt iFunc = -1) const;

			/**
			 * return the corresponding name for the input code
			 */
			void getName(const UInt& code, string& name) const;

		public:

			/**
			 * constructor for set up the information
			 */
			XCEnergyInfor(const XCIntJobInfor& infor, const XCFunc& xcfunc);

			/**
			 * destructor
			 */
			~XCEnergyInfor() { };

			///
			/// because different batch may have different length of grid points,
			/// you need to do initialization before any data collecting work!
			///
			void init(const BatchGrid& batchGrid);

			///
			/// adding results together, used for multi-threading module
			///
			void add(const XCEnergyInfor& infor) {
				for(UInt i=0; i<xcEList.size(); i++) {
					xcEList[i] += infor.xcEList[i];
				}
			};

			/**
			 * do we do xc energy profiling?
			 */
			bool doXCEnergyProfiling() const { return (nEComponents > 0); };

			///
			/// from the input weights for batch, form the batch energy and add it to the batch result
			///
			void formBatchXCEnergy(const BatchGrid& batchGrid);

			/**
			 * update the batch energy according to it's name
			 */
			void updateEComponent(const Double* F, const string& funcName, UInt iFunc = -1);

			/**
			 * print out the energy result
			 */
			void print() const; 
	};

}

#endif
