/**
 * \file    wfn.h
 * \brief   generating the wfn format of file for AIMPAC program 
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef WFN_H
#define WFN_H
#include "libgen.h"

namespace shell {
	class MolShell;
}

namespace molecule {
	class Molecule;
}

namespace mo {
	class MO;
}

namespace globalinfor {
	class GlobalInfor;
}

namespace wfn {

	using namespace mo;
	using namespace shell;
	using namespace molecule;
	using namespace globalinfor; 

	/**
	 * \class   WFN
	 * \brief   generate wfn format file
	 */
	class WFN {

		private:

			string outputFile;        ///< the name and location of output files
			UInt section;             ///< geometry section number

		public:

			/**
			 * read in the global setting
			 */
			WFN(const GlobalInfor& infor, const Molecule& mol);

			/**
			 * destructor
			 */
			~WFN() { };

			///
			/// do we really dump the wfn file?
			///
			bool doWFNFile() const {
				if (section == static_cast<UInt>(-1)) return false;
				return true;
			};

			///
			/// dump the wfn file
			///
			void dump(const Double& totalEnergy, const GlobalInfor& infor, const Molecule& mol, 
					const MolShell& ms, const MO& mo) const;
	};

}

#endif
