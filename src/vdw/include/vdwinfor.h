/**
 * \file    vdwinfor.h
 * \brief   top class to mornitor vdw calculation
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef VDW_H
#define VDW_H
#include "libgen.h"
#include <string>

namespace molecule {
	class Molecule;
}

namespace denmtrx {
	class DenMtrx;
}

namespace shell {
	class MolShell;
}

namespace xcintsinfor {
	class XCIntJobInfor;
}

namespace vdwinfor {

	using namespace molecule;
	using namespace shell;
	using namespace denmtrx;
	using namespace xcintsinfor;
	using namespace std;

	class VDWInfor {

		private:

			//
			// basic setting of VDW
			//
			UInt order;                 ///< do we do C6, C8, C10 calculations?
			bool withNonAdditiveB3;     ///< whether we add in non-additive 3 body VDW term?
			Double energy;              ///< vdw energy

			//
			// XDM section
			//
			Double a1;                  ///< first parameter - this is a radio and unitless
			Double a2;                  ///< second parameter - this is in unit of angstrom                 
			string funcName;            ///< the exchange hole functional name used in XDM module

			//
			// debug section
			//
			bool printXDMTiming;        ///< print out the timing of XDM part

		public:

			/**
			 * constructor
			 */
			VDWInfor(const GlobalInfor& input, const Molecule& mol);

			///
			/// destrucotr
			///
			~VDWInfor() { };

			///
			/// perform XDM model to get the energy
			///
			void doXDM(const Molecule& mol, const MolShell& ms, 
					const XCIntJobInfor& xcinfor, const DenMtrx& den);

			///
			/// get the energy
			///
			Double vdwEnergy() const { return energy; };

			///
			/// get the VDW order
			///
			UInt getVDWOrder() const { return order; };

			///
			/// get the XDM parameter a1
			///
			Double getXDMPar1() const { return a1; };

			///
			/// get the XDM parameter a2
			///
			Double getXDMPar2() const { return a2; };

			///
			/// whether we have non-addative B3 for XDM?
			///
			bool hasNonAddative3BodyTerm() const { return withNonAdditiveB3; };

			///
			/// get the functional name for the exchange hole for XDM
			///
			string getEXHoleFuncName() const { return funcName; };
	};
}

#endif

