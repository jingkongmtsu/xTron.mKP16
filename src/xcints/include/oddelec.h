/**
 * \file    oddelec.h
 * \brief   performing the ODD electron calculation based on the hyper functional methods
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef ODDELEC_H
#define ODDELEC_H
#include "libgen.h"
#include "scalablevec.h"

namespace shell {
	class MolShell;
}

namespace xcvar {
	class XCVar;
}

namespace batchvar {
	class BatchVar;
}

namespace denmtrx {
	class DenMtrx;
}

namespace batchgrid {
	class BatchGrid;
}

namespace xcintsinfor {
	class XCIntsinfor;
}

namespace oddelec {

	using namespace shell;
	using namespace xcvar;
	using namespace denmtrx;
	using namespace batchvar;
	using namespace batchgrid;
	using namespace xcintsinfor;

	/**
	 * \class   ODDElec
	 * \brief   perform odd electron calculation
	 *
	 * see the paper of:
	 *
	 * "Analyzing effects of strong electron correlation within Kohn-Sham density-functional theory"
	 * Emil Proynov, Fenglai Liu, and Jing Kong
	 * Phys. Rev. A 88, 032510
	 */
	class ODDElec {

		private:

			DoubleVec oddElecSum;     ///< odd electron atomic population

		public:

			/**
			 * constructor to initialize the odd electron data
			 */
			ODDElec(const MolShell& ms, const XCIntJobInfor& infor);

			/**
			 * destructor
			 */
			~ODDElec() { };

			/**
			 * perform the odd electron calculation, and sum up the result and assign it 
			 * to the given center
			 */
			void doOddElecPopulation(const XCIntJobInfor& xcinfor, const BatchGrid& bg, const BatchVar& bvar, 
					const XCVar& xcvar, const DenMtrx& den); 

			///
			/// merge results of this object with another one
			///
			void add(const XCIntJobInfor& infor, const ODDElec& oddElec);

			/**
			 * print the results
			 */
			void printResults(const MolShell& ms) const;

	};

}

#endif
